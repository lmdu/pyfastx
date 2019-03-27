#
# Cython declaration for the zran library.
#
# Author: Paul McCarthy <pauldmccarthy@gmail.com>
#

from libc.stdio  cimport FILE
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int64_t
from posix.types cimport off_t


cdef extern from "zran.h":

    ctypedef struct zran_index_t:
        FILE         *fd;
        size_t        compressed_size;
        size_t        uncompressed_size;
        uint32_t      spacing;
        uint32_t      window_size;
        uint32_t      npoints;
        zran_point_t *list;

    ctypedef struct zran_point_t:
        uint64_t  cmp_offset;
        uint64_t  uncmp_offset;
        uint8_t   bits;
        uint8_t  *data;

    enum:
        ZRAN_AUTO_BUILD       =  1,
        ZRAN_SEEK_FAIL        = -1,
        ZRAN_SEEK_OK          =  0,
        ZRAN_SEEK_NOT_COVERED =  1,
        ZRAN_SEEK_EOF         =  2,

        ZRAN_READ_NOT_COVERED = -1,
        ZRAN_READ_EOF         = -2,
        ZRAN_READ_FAIL        = -3,

        ZRAN_EXPORT_OK          =  0,
        ZRAN_EXPORT_WRITE_ERROR = -1,

        ZRAN_IMPORT_OK             =  0,
        ZRAN_IMPORT_FAIL           = -1,
        ZRAN_IMPORT_EOF            = -2,
        ZRAN_IMPORT_READ_ERROR     = -3,
        ZRAN_IMPORT_OVERFLOW       = -4,
        ZRAN_IMPORT_INCONSISTENT   = -5,
        ZRAN_IMPORT_MEMORY_ERROR   = -6,
        ZRAN_IMPORT_UNKNOWN_FORMAT = -7

    bint zran_init(zran_index_t *index,
                   FILE         *fd,
                   uint32_t      spacing,
                   uint32_t      window_size,
                   uint32_t      readbuf_size,
                   uint16_t      flags)

    void zran_free(zran_index_t *index)

    bint zran_build_index(zran_index_t *index,
                          uint64_t      from_,
                          uint64_t      until) nogil;

    long zran_tell(zran_index_t *index);

    int zran_seek(zran_index_t  *index,
                  int64_t        offset,
                  uint8_t        whence,
                  zran_point_t **point) nogil;

    int64_t zran_read(zran_index_t *index,
                      void         *buf,
                      uint64_t      len) nogil;

    int zran_export_index(zran_index_t *index,
                          FILE         *fd);

    int zran_import_index(zran_index_t *index,
                          FILE         *fd);
