# cython: binding=True,embedsignature=True
#
# The IndexedGzipFile class.
#
# Author: Paul McCarthy <pauldmccarthy@gmail.com>
#

"""This module provides the :class:`IndexedGzipFile` class, a drop-in
replacement for the built-in ``gzip.GzipFile`` class, for faster read-only
random access to gzip files.
"""


from libc.stdio     cimport (SEEK_SET,
                             SEEK_CUR,
                             FILE,
                             fopen,
                             fdopen,
                             fclose)

from libc.stdint    cimport (uint8_t,
                             uint64_t,
                             int64_t)

from cpython.mem    cimport (PyMem_Malloc,
                             PyMem_Realloc,
                             PyMem_Free)

from cpython.buffer cimport (PyObject_GetBuffer,
                             PyBuffer_Release,
                             PyBUF_ANY_CONTIGUOUS,
                             PyBUF_SIMPLE)

cimport indexed_gzip.zran as zran

import io
import warnings
import contextlib
import threading
import logging


log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


class NotCoveredError(Exception):
    """Exception raised by the :class:`_IndexedGzipFile` when an attempt is
    made to seek to/read from a location that is not covered by the
    index. This exception will never be raised if the ``_IndexedGzipFile`` was
    created with ``auto_build=True``.
    """
    pass


class ZranError(Exception):
    """Exception raised by the :class:`_IndexedGzipFile` when the ``zran``
    library signals an error.
    """
    pass


class NoHandleError(Exception):
    """Exception raised by the :class:`_IndexedGzipFile` when
    ``drop_handles is True`` and an attempt is made to access the underlying
    file object.
    """


cdef class _IndexedGzipFile:
    """The ``_IndexedGzipFile`` class allows for fast random access of a gzip
    file by using the ``zran`` library to build and maintain an index of seek
    points into the file.

    .. note:: The :meth:`seek` and :meth:`read` methods release the GIL while
              calling ``zran`` functions, but the ``_IndexedGzipFile`` is *not*
              thread-safe. Use the ``IndexedGzipFile`` class (i.e. without the
              leading underscore) if you need thread-safety.
    """


    cdef zran.zran_index_t index
    """A reference to the ``zran_index`` struct. """


    cdef readonly bint auto_build
    """Flag which is set to ``True`` if the file index is built automatically
    on seeks/reads.
    """


    cdef readonly unsigned int readall_buf_size
    """Used by meth:`read` as the read buffer size."""


    cdef readonly object filename
    """String containing path of file being indexed. Used to release and
    reopen file handles between seeks and reads.
    Set to ``None`` if file handle is passed.
    """


    cdef readonly bint own_file
    """Flag which tracks whether this ``_IndexedGzipFile`` has opened its
    own file handle, or was given one.
    """


    cdef readonly bint drop_handles
    """Copy of the ``drop_handles`` flag as passed to :meth:`__cinit__`. """


    cdef object pyfid
    """A reference to the python file handle. """


    cdef bint finalized
    """Flag which is set to ``True`` if the ``_IndexedGzipFile`` has been
    closed. Further operations will fail if ``True``.
    """


    def __cinit__(self,
                  filename=None,
                  fid=None,
                  fileobj=None,
                  mode=None,
                  auto_build=True,
                  spacing=4194304,
                  window_size=32768,
                  readbuf_size=1048576,
                  readall_buf_size=16777216,
                  drop_handles=True,
                  index_file=None):
        """Create an ``_IndexedGzipFile``. The file may be specified either
        with an open file handle (``fileobj``), or with a ``filename``. If the
        former, the file must have been opened in ``'rb'`` mode.

        :arg filename:         File name.

        :arg mode:             Opening mode. Must be either ``'r'`` or ``'rb``.

        :arg fileobj:          Open file handle.

        :arg fid:              Deprecated, use ``fileobj`` instead.

        :arg auto_build:       If ``True`` (the default), the index is
                               automatically built on seeks/reads.

        :arg spacing:          Number of bytes between index seek points.

        :arg window_size:      Number of bytes of uncompressed data stored with
                               each seek point.

        :arg readbuf_size:     Size of buffer in bytes for storing compressed
                               data read in from the file.

        :arg readall_buf_size: Size of buffer in bytes used by :meth:`read`
                               when reading until EOF.

        :arg drop_handles:     Has no effect if an open ``fid`` is specified,
                               rather than a ``filename``.  If ``True`` (the
                               default), a handle to the file is opened and
                               closed on every access. Otherwise the file is
                               opened at ``__cinit__``, and kept open until
                               this ``_IndexedGzipFile`` is destroyed.

        :arg index_file:       Pre-generated index for this ``gz`` file -
                               if provided, passed through to
                               :meth:`import_index`.
        """

        cdef FILE *fd = NULL

        if fid is not None:
            warnings.warn('fid is deprecated - use fileobj instead',
                          DeprecationWarning)

        if fid is not None and fileobj is None:
            fileobj = fid
            fid     = None

        if (fileobj is     None and filename is     None) or \
           (fileobj is not None and fid      is not None) or \
           (fileobj is not None and filename is not None):
            raise ValueError('One of fileobj or filename must be specified')

        if fileobj is not None and fileobj.mode != 'rb':
            raise ValueError('The gzip file must be opened in '
                             'read-only binary ("rb") mode')

        if (fileobj is None) and (mode not in (None, 'r', 'rb')):
            raise ValueError('Invalid mode ({}), must be '
                             '\'r\' or \'rb\''.format(mode))

        mode     = 'rb'
        own_file = fileobj is None

        # if file is specified with an open
        # file handle, drop_handles is ignored
        if fileobj is not None:
            drop_handles = False

        # if not drop_handles, we open a
        # file handle and keep it open for
        # the lifetime of this object.
        if not drop_handles:
            if fileobj is None:
                fileobj = open(filename, mode)
            fd = fdopen(fileobj.fileno(), 'rb')

        self.auto_build       = auto_build
        self.readall_buf_size = readall_buf_size
        self.drop_handles     = drop_handles
        self.filename         = filename
        self.own_file         = own_file
        self.pyfid            = fileobj
        self.index.fd         = fd

        if self.auto_build: flags = zran.ZRAN_AUTO_BUILD
        else:               flags = 0

        with self.__file_handle():
            if zran.zran_init(index=&self.index,
                              fd=self.index.fd,
                              spacing=spacing,
                              window_size=window_size,
                              readbuf_size=readbuf_size,
                              flags=flags):
                raise ZranError('zran_init returned error')

        log.debug('{}.__init__({}, {}, {}, {}, {}, {}, {})'.format(
            type(self).__name__,
            fileobj,
            filename,
            auto_build,
            spacing,
            window_size,
            readbuf_size,
            drop_handles))

        if index_file is not None:
            self.import_index(index_file)


    def __init__(self, *args, **kwargs):
        """This method does nothing. It is here to make sub-classing
        ``_IndexedGzipFile`` easier.
        """
        pass


    def __file_handle(self):
        """This method is used as a context manager whenever access to the
        underlying file stream is required. It makes sure that ``index.fd``
        field is set appropriately, opening/closing the file handle as
        necessary (depending on the value of :attr:`drop_handles`).
        """

        # Errors occur with Python 2.7 and
        # Cython < 0.26 when decorating
        # cdef-class methods. This workaround
        # can be removed when you are happy
        # dropping support for cython < 0.26.
        @contextlib.contextmanager
        def proxy():

            # If a file handle already exists,
            # return it. This clause makes this
            # context manager reentrant.
            if self.index.fd is not NULL:
                yield

            # otherwise we open a new
            # file handle on each access
            else:

                try:
                    self.index.fd = fopen(self.filename.encode(), 'rb')
                    yield

                finally:
                    fclose(self.index.fd)
                    self.index.fd = NULL

        return proxy()


    def fileno(self):
        """Calls ``fileno`` on the underlying file object. Raises a
        :exc:`NoHandleError` if ``drop_handles is True``.
        """
        if self.drop_handles:
            raise NoHandleError()
        return self.pyfid.fileno()


    def fileobj(self):
        """Returns a reference to the python file object. Raises a
        :exc:`NoHandleError` if ``drop_handles is True``.
        """
        if self.drop_handles:
            raise NoHandleError()
        return self.pyfid


    @property
    def mode(self):
        """Returns the mode that this file was opened in. Currently always
        returns ``'rb'``.
        """
        return 'rb'


    def close(self):
        """Closes this ``_IndexedGzipFile``. """

        if self.closed:
            raise IOError('_IndexedGzipFile is already closed')

        if   self.own_file and self.pyfid    is not None: self.pyfid.close()
        elif self.own_file and self.index.fd is not NULL: fclose(self.index.fd)

        zran.zran_free(&self.index)

        self.index.fd  = NULL
        self.filename  = None
        self.pyfid     = None
        self.finalized = True

        log.debug('{}.close()'.format(type(self).__name__))


    @property
    def closed(self):
        """Returns ``True`` if this ``_IndexedGzipFile`` is closed, ``False``
        otherwise.
        """
        return self.finalized


    def readable(self):
        """Returns ``True`` if this ``_IndexedGzipFile`` is readable, ``False``
        otherwise.
        """
        return not self.closed


    def writable(self):
        """Currently always returns ``False`` - the ``_IndexedGzipFile`` does
        not support writing yet.
        """
        return False


    def seekable(self):
        """Returns ``True`` if this ``_IndexedGzipFile`` supports seeking,
        ``False`` otherwise.
        """
        return not self.closed


    def tell(self):
        """Returns the current seek offset into the uncompressed data stream.
        """
        return zran.zran_tell(&self.index)


    def __enter__(self):
        """Returns this ``_IndexedGzipFile``. """
        return self


    def __exit__(self, *args):
        """Calls close on this ``_IndexedGzipFile``. """
        if not self.closed:
            self.close()


    def __dealloc__(self):
        """Frees the memory used by this ``_IndexedGzipFile``. If a file name
        was passed to :meth:`__cinit__`, the file handle is closed.
        """
        if not self.closed:
            self.close()


    def build_full_index(self):
        """Re-builds the full file index. """

        with self.__file_handle():
            ret = zran.zran_build_index(&self.index, 0, 0)

        if ret != 0:
            raise ZranError('zran_build_index returned error')

        log.debug('{}.build_full_index()'.format(type(self).__name__))


    def seek(self, offset, whence=SEEK_SET):
        """Seeks to the specified position in the uncompressed data stream.

        If this ``_IndexedGzipFile`` was created with ``auto_build=False``,
        and the requested offset is not covered by the index, a
        :exc:`NotCoveredError` is raised.

        :arg offset: Desired seek offset into the uncompressed data

        :arg whence: Either  ``SEEK_SET`` or ``SEEK_CUR``. If not one of these,
                     a :exc:`ValueError` is raised.

        :returns:    The final seek location into the uncompressed stream.

        .. note:: This method releases the GIL while ``zran_seek`` is
                  running.
        """

        cdef int                ret
        cdef int64_t            off      = offset
        cdef uint8_t            c_whence = whence
        cdef zran.zran_index_t *index    = &self.index

        if whence not in (SEEK_SET, SEEK_CUR):
            raise ValueError('Seek from end not supported')

        with self.__file_handle(), nogil:
            ret = zran.zran_seek(index, off, c_whence, NULL)

        if ret < 0:
            raise ZranError('zran_seek returned error: {}'.format(ret))

        elif ret == zran.ZRAN_SEEK_NOT_COVERED:
            raise NotCoveredError('Index does not cover '
                                  'offset {}'.format(offset))

        elif ret not in (zran.ZRAN_SEEK_OK, zran.ZRAN_SEEK_EOF):
            raise ZranError('zran_seek returned unknown code: {}'.format(ret))

        offset = self.tell()

        log.debug('{}.seek({})'.format(type(self).__name__, offset))

        return offset


    def read(self, nbytes=-1):
        """Reads up to ``nbytes`` bytes from the uncompressed data stream.
        If ``nbytes < 0`` the stream is read until EOF.

        If the stream is already at EOF, ``b''`` is returned.

        .. note:: This method releases the GIL while ``zran_read`` is
                  running.
        """

        if   nbytes == 0: return bytes()
        elif nbytes <  0: buf = ReadBuffer(self.readall_buf_size)
        else:             buf = ReadBuffer(nbytes)

        cdef zran.zran_index_t *index  = &self.index
        cdef size_t             nread  = 0
        cdef uint64_t           bufsz  = buf.size
        cdef size_t             offset = 0
        cdef void              *buffer
        cdef int64_t            ret

        with self.__file_handle():

            # Read until EOF or enough
            # bytes have been read
            while True:

                # read some bytes into the correct
                # buffer location
                buffer = <char *>buf.buffer + offset
                with nogil:
                    ret = zran.zran_read(index, buffer, bufsz)

                # see how the read went
                if ret == zran.ZRAN_READ_FAIL:
                    raise ZranError('zran_read returned error ({})'.format(ret))

                # This will happen if the current
                # seek point is not covered by the
                # index, and auto-build is disabled
                elif ret == zran.ZRAN_READ_NOT_COVERED:
                    raise NotCoveredError('Index does not cover current offset')

                # No bytes were read, and there are
                # no more bytes to read. This will
                # happen when the seek point was at
                # or beyond EOF when zran_read was
                # called
                elif ret == zran.ZRAN_READ_EOF:
                    break

                nread  += ret
                offset += ret

                # If we requested a specific number of
                # bytes, zran_read will have returned
                # them all (or all until EOF), so we're
                # finished
                if nbytes > 0:
                    break

                # Otherwise if reading until EOF, check
                # and increase the buffer size if necessary
                if (nread + self.readall_buf_size) > buf.size:
                    buf.resize(buf.size + self.readall_buf_size)
                    offset = nread

        buf.resize(nread)
        pybuf = <bytes>(<char *>buf.buffer)[:nread]

        log.debug('{}.read({})'.format(type(self).__name__, len(pybuf)))

        return pybuf


    def readinto(self, buf):
        """Reads up to ``len(buf)`` bytes directly into ``buf``, which is
        assumed to be a mutable ``bytes``-like object (e.g. a ``memoryview``
        or ``bytearray``.
        """

        cdef zran.zran_index_t *index  = &self.index
        cdef uint64_t           bufsz  = len(buf)
        cdef Py_buffer          pbuf
        cdef void              *vbuf
        cdef int64_t            ret

        # Create a Py_Buffer which allows
        # us to access the memory managed
        # by the provided buf
        PyObject_GetBuffer(buf, &pbuf, PyBUF_SIMPLE | PyBUF_ANY_CONTIGUOUS)

        # read some bytes
        try:

            vbuf = <void *>pbuf.buf
            with self.__file_handle(), nogil:
                ret = zran.zran_read(index, vbuf, bufsz)

        # release the py_buffer
        finally:
            PyBuffer_Release(&pbuf)

        # see how the read went
        if ret == zran.ZRAN_READ_FAIL:
            raise ZranError('zran_read returned error ({})'.format(ret))

        # This will happen if the current
        # seek point is not covered by the
        # index, and auto-build is disabled
        elif ret == zran.ZRAN_READ_NOT_COVERED:
            raise NotCoveredError('Index does not cover current offset')

        # No bytes were read, and there are
        # no more bytes to read. This will
        # happen when the seek point was at
        # or beyond EOF when zran_read was
        # called
        elif ret == zran.ZRAN_READ_EOF:
            return 0

        # Return the number of bytes that
        # were read
        else:
            return ret


    def pread(self, nbytes, offset):
        """Seeks to the specified ``offset``, then reads and returns
        ``nbytes``. See :meth:`seek` and :meth:`read`.
        """
        with self.__file_handle():
            self.seek(offset)
            return self.read(nbytes)


    def readline(self, size=-1):
        """Read and return up to the next ``'\n'`` character (up to at most
        ``size`` bytes, if ``size >= 0``) from the uncompressed data stream.

        If the end of the stream has been reached, ``b''`` is returned.
        """

        if size == 0:
            return bytes()

        linebuf  = b''
        startpos = self.tell()
        bufsz    = 1024

        # Read in chunks of [bufsz] bytes at a time
        with self.__file_handle():
            while True:

                buf = self.read(bufsz)

                lineidx  = buf.find(b'\n')
                haveline = lineidx  >= 0
                eof      = len(buf) == 0

                # Are we at EOF? Nothing more to do
                if eof:
                    break

                # Have we found a line? Discard
                # everything that comes after it
                if haveline:
                    linebuf = linebuf + buf[:lineidx + 1]

                # If we've found a line, and are
                # not size-limiting, we're done
                if haveline and size < 0:
                    break

                # If we're size limiting, and have
                # read in enough bytes, we're done
                if size >= 0 and len(linebuf) > size:
                    linebuf = linebuf[:size]
                    break

            # Rewind the seek location
            # to the finishing point
            self.seek(startpos + len(linebuf))

        return linebuf


    def readlines(self, hint=-1):
        """Reads and returns a list of lines from the uncompressed data.
        If ``hint`` is provided, lines will be read until the total size
        of all lines exceeds ``hint`` in bytes.
        """

        totalsize = 0
        lines     = []

        with self.__file_handle():
            while True:

                line = self.readline()
                if line == b'':
                    break

                lines.append(line)

                totalsize += len(line)

                if hint >= 0 and totalsize > hint:
                    break

        return lines


    def __iter__(self):
        """Returns this ``_IndexedGzipFile`` which can be iterated over to
        return lines (separated by ``'\n'``) in the uncompressed stream.
        """
        return self


    def __next__(self):
        """Returns the next line from the uncompressed stream. Raises
        :exc:`StopIteration` when there are no lines left.
        """
        line = self.readline()

        if line == b'':
            raise StopIteration()
        else:
            return line


    def write(self, *args, **kwargs):
        """Currently raises a :exc:`NotImplementedError`."""
        raise NotImplementedError('_IndexedGzipFile does not support writing')


    def flush(self):
        """Currently does nothing. """
        pass


    def export_index(self, filename=None, fileobj=None):
        """Export index data to the given file. Either ``filename`` or
        ``fileobj`` should be specified, but not both. ``fileobj`` should be
        opened in 'wb' mode.

        :arg filename: Name of the file.
        :arg fileobj:  Open file handle.
        """

        if filename is None and fileobj is None:
            raise ValueError('One of filename or fileobj must be specified')

        if filename is not None and fileobj is not None:
            raise ValueError(
                'Only one of filename or fileobj must be specified')

        if filename is not None:
            fileobj = open(filename, 'wb')
            close_file = True

        else:
            close_file = False
            if fileobj.mode != 'wb':
                raise ValueError(
                    'File should be opened write-only binary mode.')

        try:
            fd  = fdopen(fileobj.fileno(), 'ab')
            ret = zran.zran_export_index(&self.index, fd)
            if ret != zran.ZRAN_EXPORT_OK:
                raise ZranError('export_index returned error: {}', ret)

        finally:
            if close_file:
                fileobj.close()

        log.debug('{}.export_index({}, {})'.format(
            type(self).__name__,
            filename,
            fileobj))


    def import_index(self, filename=None, fileobj=None):
        """Import index data from the given file. Either ``filename`` or
        ``fileobj`` should be specified, but not both. ``fileobj`` should be
        opened in 'rb' mode.

        :arg filename: Name of the file.
        :arg fileobj:  Open file handle.
        """

        if filename is None and fileobj is None:
            raise ValueError('One of filename or fileobj must be specified')

        if filename is not None and fileobj is not None:
            raise ValueError(
                'Only one of filename or fileobj must be specified')

        if filename is not None:
            fileobj    = open(filename, 'rb')
            close_file = True

        else:
            close_file = False
            if fileobj.mode != 'rb':
                raise ValueError(
                    'File should be opened read-only binary mode.')

        try:
            fd  = fdopen(fileobj.fileno(), 'rb')
            ret = zran.zran_import_index(&self.index, fd)
            if ret != zran.ZRAN_IMPORT_OK:
                raise ZranError('import_index returned error: {}', ret)

        finally:
            if close_file:
                fileobj.close()

        log.debug('{}.import_index({}, {})'.format(
            type(self).__name__,
            filename,
            fileobj))


cdef class ReadBuffer:
    """Wrapper around a chunk of memory.

    .. see:: http://docs.cython.org/src/tutorial/memory_allocation.html
    """

    cdef void *buffer
    """A raw chunk of bytes. """


    cdef size_t size;
    """Size of the buffer. """


    def __cinit__(self, size_t size):
        """Allocate ``size`` bytes of memory. """

        self.size   = size
        self.buffer = PyMem_Malloc(size)

        if not self.buffer:
            raise MemoryError('PyMem_Malloc fail')

        log.debug('ReadBuffer.__cinit__({})'.format(size))


    def resize(self, size_t size):
        """Re-allocate the memory to the given ``size``. """

        if size == self.size:
            return

        buf = PyMem_Realloc(self.buffer, size)

        if not buf:
            raise MemoryError('PyMem_Realloc fail')

        log.debug('ReadBuffer.resize({})'.format(size))

        self.size   = size
        self.buffer = buf


    def __dealloc__(self):
        """Free the mwmory. """
        PyMem_Free(self.buffer)

        log.debug('ReadBuffer.__dealloc__()')


class IndexedGzipFile(io.BufferedReader):
    """The ``IndexedGzipFile`` is an ``io.BufferedReader`` which wraps
    an :class:`_IndexedGzipFile` instance. By accessing the
    ``_IndexedGzipFile`` instance through an ``io.BufferedReader``, read
    performance is improved through buffering, and access to the I/O methods
    is made thread-safe.

    A :meth:`pread` method is also implemented, as it is not implemented by
    the ``io.BufferedReader``.
    """


    def __init__(self, *args, **kwargs):
        """Opens an ``_IndexedGzipFile``, and then calls
        ``io.BufferedReader.__init__``.

        :arg buffer_size: Optional, must be passed as a keyword argument.
                          Passed through to ``io.BufferedReader.__init__``.
                          If not provided, a default value of 1048576 is used.

        All other arguments are passed through to
        ``_IndezedGzipFile.__init__``.
        """

        buffer_size     = kwargs.pop('buffer_size', 1048576)
        fobj            = _IndexedGzipFile(*args, **kwargs)
        self.__fileLock = threading.RLock()

        self.build_full_index = fobj.build_full_index
        self.import_index     = fobj.import_index
        self.export_index     = fobj.export_index
        self.fileobj          = fobj.fileobj
        self.drop_handles     = fobj.drop_handles

        super(IndexedGzipFile, self).__init__(fobj, buffer_size)


    def pread(self, nbytes, offset):
        """Seeks to ``offset``, then reads and returns up to ``nbytes``.
        The calls to seek and read are protected by a ``threading.RLock``.
        """
        with self.__fileLock:
            self.seek(offset)
            return self.read(nbytes)


class SafeIndexedGzipFile(IndexedGzipFile):
    """Deprecated - use :class:`IndexedGzipFile` instead. """
    def __init__(self, *args, **kwargs):
        """Deprecated - use :class:`IndexedGzipFile` instead. """
        warnings.warn('SafeIndexedGzilFile is deprecated - '
                      'use IndexedGzipFile instead',
                      DeprecationWarning)
        super(SafeIndexedGzipFile, self).__init__(*args, **kwargs)
