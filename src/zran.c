/*
 * zran.c - indexed access to gzip files.
 *
 * See zran.h for documentation.
 *
 * This module was originally based on the zran example, written by Mark
 * Adler, which ships with the zlib source code.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#ifdef _WIN32
#include "windows.h"
#include "io.h"
static int is_readonly(FILE *fd, PyObject *f)
{
    /* Can't find a way to do this correctly under
       Windows and the check is not required anyway
       since the underlying Python module checks it
       already */
    return 1;
}
#else
#include <fcntl.h>
/* Check if file is read-only */
static int is_readonly(FILE *fd, PyObject *f)
{
    /* Skip the test for python file-likes */
    return fd != NULL
      ? (fcntl(fileno(fd), F_GETFL) & O_ACCMODE) == O_RDONLY
      : 1;
}

static uint32_t max(uint32_t a, uint32_t b) {

  if (a > b) return a;
  else       return b;
}
#endif


#include "zran.h"
#include "zran_file_util.h"


#ifdef NO_C99
static double round(double val)
{
    return floor(val + 0.5);
}
#endif

/*
 * Turn this on to make noise.
 *
 * #define ZRAN_VERBOSE
 */
//#define ZRAN_VERBOSE


#ifdef ZRAN_VERBOSE
#define zran_log(...) fprintf(stderr, __VA_ARGS__)
#else
#define zran_log(...)
#endif


/*
 * Identifier and version number for index files created by zran_export_index.
 */
const char    ZRAN_INDEX_FILE_ID[]    = {'G', 'Z', 'I', 'D', 'X'};
const uint8_t ZRAN_INDEX_FILE_VERSION = 1;


/*
 * Discards all points in the index which come after the specified
 * compressed offset.
 *
 * Returns 0 on success, non-0 on failure.
 */
static int _zran_invalidate_index(
    zran_index_t *index, /* The index                       */
    uint64_t      from   /* Offset into the compressed data */
);


/*
 * Expands the capacity of the memory used to store the index list.
 *
 * Returns 0 on success, non-0 on failure.
 */
static int _zran_expand_point_list(
    zran_index_t *index /* The index */
);


/*
 * Reduces the capacity of the memory used to store the index list, so that it
 * is only as big as necessary.
 *
 * Returns 0 on success, non-0 on failure.
 */
static int _zran_free_unused(
    zran_index_t *index  /* The index */
);

/*
 * Returns the current limit of the index, i.e. how much of the file is covered
 * by the index.
 */
static uint64_t _zran_index_limit(
  zran_index_t *index,      /* The index */
  uint8_t       compressed  /* Pass in non-0 to get the compressed stream
                               limit, or 0 for the uncompressed limit. */
);


/* Return codes for _zran_get_point_at */
int ZRAN_GET_POINT_CRC_ERROR   =  -2;
int ZRAN_GET_POINT_FAIL        =  -1;
int ZRAN_GET_POINT_OK          =   0;
int ZRAN_GET_POINT_NOT_COVERED =   1;
int ZRAN_GET_POINT_EOF         =   2;

/*
 * Searches for the zran_point which precedes the given offset. The offset
 * may be specified as being relative to the start of the compressed data,
 * or the uncompressed data.
 *
 * Returns:
 *
 *   - ZRAN_GET_POINT_OK on success.
 *
 *   - ZRAN_GET_POINT_NOT_COVERED if the index does not yet cover the
 *     specified offset.

 *   - ZRAN_GET_POINT_EOF if the specified offset is at or beyond the end of
 *     the file.
 */
static int _zran_get_point_at(
    zran_index_t  *index,      /* The index */

    uint64_t       offset,     /* The desired offset into the compressed or
                                  uncompressed data stream */

    uint8_t        compressed, /* Pass in 0 or non-0 to indicate that the
                                  offset is relative to the uncompressed or
                                  compressed data streams, respectively. */

    zran_point_t **point       /* If an index point corresponding to the
                                  specified offset is identified, this pointer
                                  will be updated to point to it. */
);


/*
 * If the index has been created without the ZRAN_AUTO_BUILD flag, this
 * function is identical to the _zran_get_point_at function.
 *
 * If the index has been created with the ZRAN_AUTO_BUILD flag, and the
 * requested offset is beyond the current range of the index, the index will
 * be expanded to encompass it.
 *
 * The input arguments and return values are identical to the
 * _zran_get_point_at function, however:
 *
 *   - if the index has been initialised with the ZRAN_AUTO_BUILD flag, this
 *     function will never return ZRAN_GET_POINT_NOT_COVERED.
 *
 *   - If a CRC validation error occurs while the index is being expanded,
 *     ZRAN_GET_POINT_CRC_ERROR is returned.
 */
static int _zran_get_point_with_expand(
    zran_index_t  *index,      /* The index                           */
    uint64_t       offset,     /* Desired offset                      */
    uint8_t        compressed, /* Compressed or uncompressed offset   */
    zran_point_t **point       /* Place to store the identified point */
);


/*
 * Estimate an offset in the compressed / uncompressed data stream
 * corresponding to the given offset, which is specified in the uncompressed /
 * compressed data stream.  If the given offset is specified relative to the
 * compressed data stream, the returned value is a location in the
 * uncompressed data stream which approximately corresponds to the given
 * offset.
 *
 * This function is used by the _zran_get_point_with_expand function, if the
 * index has been created with the ZRAN_AUTO_BUILD flag, to determine how far
 * the index needs to be expanded to cover a requested offset that is not yet
 * covered.
 */
static uint64_t _zran_estimate_offset(
    zran_index_t *index,      /* The index */

    uint64_t      offset,     /* The offset for which a corresponding offset
                                 is to be estimated. */

    uint8_t       compressed  /* Pass in 0 or non-0 to indicate that the given
                                 offset is specified relative to the
                                 uncompressed or compressed stream,
                                 respectively. */
);


/*
 * Used by _zran_inflate. Initialises zlib to start decompressing/inflating
 * from either:
 *
 *  - the current seek location in the compressed data, or
 *  - from a location denoted by a specific index point
 *
 * If an index point is provided, the function will seek to the specified
 * compressed data offset before initialising zlib.
 *
 * Otherwise (no index point), inflation is initialised at the current seek
 * location in the input data, and a GZIP header is expected at that location.
 *
 * The index->readbuf and readbuf_size, and the z_stream->avail_in, avail_out,
 * next_in and next_out fields must all be set before this function is called.
 *
 * Returns the number of bytes over the input data that were read (which could
 * be 0), or a negative value on failure.
 */
static int _zran_init_zlib_inflate(
    zran_index_t *index,        /* The index */

    z_stream     *stream,       /* Pointer to a z_stream struct */

    zran_point_t *point         /* Pass in NULL to initialise for inflation
                                   from the current location in the input file.
                                   Or pass a pointer to the index point
                                   corresponding to the location to start
                                   from. */
);


/*
 * Return codes for _zran_expand_index. These are currently
 * assumed to have identical values to the ZRAN_BUILD_INDEX
 * return codes.
 */
int ZRAN_EXPAND_INDEX_OK        =  0;
int ZRAN_EXPAND_INDEX_FAIL      = -1;
int ZRAN_EXPAND_INDEX_CRC_ERROR = -2;


/*
 * Expands the index from its current end-point until the given offset (which
 * must be specified relative to the compressed data stream).
 *
 * The index is expanded so that the last point comes after the given offset.
 * If the specified offset is past the last point in the index, a call to this
 * function is guaranteed to create at least one more index point. If there
 * is already an index point which comes after the offset, this function does
 * nothing, and return a success code.
 *
 * We require at least one point to be created, because we want index points
 * to be located at compression block boundaries, but in some data there may
 * be a long distance between block boundaries (longer than the desired index
 * point spacing).
 *
 * Returns 0 on success. If a CRC check fails, returns
 * ZRAN_EXPAND_INDEX_CRC_ERROR. For other types of failure, returns
 * ZRAN_EXPAND_INDEX_FAIL.
 */
static int _zran_expand_index(
    zran_index_t *index, /* The index */
    uint64_t      until  /* Expand the index to this point. If 0,
                            expand the index until EOF is reached. */
);


/*
 * Adds a new point to the end of the index.
 */
static int _zran_add_point(
    zran_index_t *index,        /* The index */

    uint8_t       bits,         /* If the compressed and uncompressed offsets
                                   are not byte-aligned, this is the number
                                   of bits in the compressed data, before the
                                   cmp_offset, where the point is located. */

    uint64_t      cmp_offset,   /* Offset into the compressed data. */

    uint64_t      uncmp_offset, /* Offset into the uncompressed data. */

    uint32_t      data_offset,  /* Offset into the data pointer specifying the
                                   point at which the uncompressed data
                                   associated with this point begins - see
                                   _zran_expand_index. It is assumed that the
                                   uncompressed data wraps around this
                                   offset. */

    uint32_t      data_size,    /* Number of bytes in data */

    uint8_t      *data          /* Pointer to data_size bytes of uncompressed
                                   data preceding this index point. */
);


/* _zran_read_data return codes */
int ZRAN_READ_DATA_EOF   = -1;
int ZRAN_READ_DATA_ERROR = -2;

/*
 * This function is a sub-function of _zran_inflate, used to read data from
 * the input file to be passed to zlib:inflate for decompression.
 *
 * Up to index->readbuf_size bytes are read from the input file into
 * index->readbuf, and the z_stream counters/pointers are updated accordingly.
 *
 * On success, returns 0.
 *
 * If there are no more bytes left to read from the input file (i.e. we are at
 * EOF), ZRAN_READ_DATA_EOF is returned. If an error occurs, returns
 * ZRAN_READ_DATA_ERROR.
 */
static int _zran_read_data_from_file(
    zran_index_t *index,        /* The index                               */
    z_stream     *stream,       /* The z_stream struct                     */
    uint64_t      cmp_offset,   /* Current offset in the compressed data   */
    uint64_t      uncmp_offset, /* Current offset in the uncompressed data */
    uint32_t      need_atleast  /* Skip read if the read buffer already has
                                   this many bytes */
);


/* _zran_find_next_stream return codes */
int ZRAN_FIND_STREAM_ERROR     = -2;
int ZRAN_FIND_STREAM_NOT_FOUND = -1;


/*
 * This function is a sub-function of _zran_inflate, used to search for a new
 * GZIP stream in a series of concatenated streams.  It searches through the
 * compressed data (pointed to by stream->next_in) to find the location of the
 * next compressed stream.
 *
 * If a new stream was found, the z_stream struct is re-initialised to
 * decompress data from the new stream, using _zran_init_zlib_inflate. In
 * this case, the value returned by that function is returned.
 *
 * Otherwise (if a compressed stream was not found), this function returns
 * ZRAN_FIND_STREAM_NOT_FOUND.
 *
 * The number of bytes that were skipped over before the new stream was found
 * is added to the provided offset pointer.
 *
 * If an error occurs, ZRAN_FIND_STREAM_ERROR is returned.
 */
static int _zran_find_next_stream(
    zran_index_t *index,  /* The index                                      */
    z_stream     *stream, /* The z_stream struct                            */
    int          *offset  /* Used to store the number of bytes skipped over */
);


/* _zran_validate_stream return codes */
int ZRAN_VALIDATE_STREAM_ERROR   = -2;
int ZRAN_VALIDATE_STREAM_INVALID = -1;


/*
 * This function is a sub-function of _zran_inflate, called when the end of a
 * gzip stream is reached. It reads the CRC32 and uncompressed file size from
 * the end of the stream, and compares them to the CRC32 and size that was
 * incrementally calculated by _zran_inflate (which are stored in
 * index->stream_crc32 and index->stream_size),
 *
 * The number of bytes that were read before the new stream was found is
 * added to the provided offset pointer.
 *
 * If ZRAN_SKIP_CRC_CHECK is active, this function returns immediately without
 * doing anything.
 *
 * If an error occurs, ZRAN_VALIDATE_STREAM_ERROR is returned.
 */
static int _zran_validate_stream(
    zran_index_t *index,  /* The index                                      */
    z_stream     *stream, /* The z_stream struct                            */
    int          *offset  /* Used to store the number of bytes skipped over */
);


/* _zran_inflate return codes */
int ZRAN_INFLATE_CRC_ERROR      = -6;
int ZRAN_INFLATE_ERROR          = -5;
int ZRAN_INFLATE_NOT_COVERED    = -4;
int ZRAN_INFLATE_OUTPUT_FULL    = -3;
int ZRAN_INFLATE_BLOCK_BOUNDARY = -2;
int ZRAN_INFLATE_EOF            = -1;
int ZRAN_INFLATE_OK             =  0;

/*
 * _zran_inflate input flags.
 * Bit position, as a power of 2
 */
uint32_t ZRAN_INFLATE_INIT_Z_STREAM         = 1;
uint32_t ZRAN_INFLATE_FREE_Z_STREAM         = 2;
uint32_t ZRAN_INFLATE_INIT_READBUF          = 4;
uint32_t ZRAN_INFLATE_FREE_READBUF          = 8;
uint32_t ZRAN_INFLATE_USE_OFFSET            = 16;
uint32_t ZRAN_INFLATE_CLEAR_READBUF_OFFSETS = 32;
uint32_t ZRAN_INFLATE_STOP_AT_BLOCK         = 64;


/* Macros used by _zran_inflate for testing flags. */
#define inflate_init_stream(  flags) ((flags & ZRAN_INFLATE_INIT_Z_STREAM) > 0)
#define inflate_free_stream(  flags) ((flags & ZRAN_INFLATE_FREE_Z_STREAM) > 0)
#define inflate_init_readbuf( flags) ((flags & ZRAN_INFLATE_INIT_READBUF)  > 0)
#define inflate_free_readbuf( flags) ((flags & ZRAN_INFLATE_FREE_READBUF)  > 0)
#define inflate_use_offset(   flags) ((flags & ZRAN_INFLATE_USE_OFFSET)    > 0)
#define inflate_stop_at_block(flags) ((flags & ZRAN_INFLATE_STOP_AT_BLOCK) > 0)
#define inflate_clear_readbuf_offsets(flags) \
    ((flags & ZRAN_INFLATE_CLEAR_READBUF_OFFSETS) > 0)


/*
 * Inflate (decompress) the specified number of bytes, or until the next
 * Z_BLOCK/Z_STREAM_END is reached.
 *
 * This is a complicated function which implements the core decompression
 * routine, and is used by both _zran_expand_index, and zran_read. It reads
 * compressed data from the file, starting from the specified compressed
 * offset, inflates (a.k.a. decompresses) it, and copies the decompressed
 * data to the provided output buffer.
 *
 * This function is complicated because it is used in three different
 * situations:
 *   - When generating the index (by zran_expand_index)
 *   - When starting from an index seek point and discarding compressed data
 *     to find a requested seek location (by zran_read)
 *   - When actually reading and decompressing data (by zran_read).
 *
 * This function may be used in a re-entrant or non-re-entrant manner,
 * depending on the flags which are used. In the latter (more likely) case,
 * various pieces of information representing the current inflation state are
 * stored in fields of the zran_index_t struct.
 *
 * Specifically, this function does the following:
 *
 *   1. Figures out the starting offsets into the compressed/uncompressed
 *      streams. If the ZRAN_INFLATE_USE_OFFSET flag is active, the index
 *      point preceding the specified offset is used as the starting point.
 *      If there is no such point, ZRAN_INFLATE_NOT_COVERED is returned.  If
 *      ZRAN_INFLATE_USE_OFFSET is not active, index->inflate_cmp_offset and
 *      index->inflate_uncmp_offset are used as the starting point.
 *
 *   2. Create a read buffer, if ZRAN_INFLATE_INIT_READBUF is active. A
 *      reference to the read buffer is stored at index->readbuf.  If
 *      ZRAN_INFLATE_INIT_READBUF is not set, the function assumes that the
 *      read buffer already exists.
 *
 *   3. If the ZRAN_INFLATE_CLEAR_READBUF_OFFSETS flag is active, the read
 *      buffer offset (index->readbuf_offset) and length (index->readbuf_end)
 *      fields are both set to 0. Otherwise, the function assumes that the
 *      current offset/length values are valid.

 *   4. Initialises the z_stream struct, if ZRAN_INFLATE_INIT_Z_STREAM is
 *      active. Otherwise, the function assumes that the z_stream struct is
 *      already initialised and ready to be used.
 *
 *   5. Read some compressed data from the file into the read buffer as needed.
 *
 *   6. Pass that data to the zlib inflate function, and store the resulting
 *      uncompressed data in the provided data buffer. If the end of a GZIP
 *      stream is reached for the first time, it is validated against the
 *      CRC/file size stored in the GZIP footer (unless ZRAN_SKIP_CRC_CHECK
 *      is active).
 *
 *   7. Repeat steps 5 and 6 until one of the following is true:
 *
 *       - The requested number of bytes have been read
 *
 *       - The output buffer is full
 *
 *       - End of file is reached
 *
 *       - ZRAN_INFLATE_STOP_AT_BLOCK is active, and a block is reached
 *
 *   8. If ZRAN_INFLATE_FREE_READBUF is active, the file read buffer is
 *      de-allocated.
 *
 *   9. If ZRAN_INFLATE_FREE_Z_STREAM is active, the memory used by the
 *      z_stream struct is de-allocated (via the zlib inflateEnd function).
 *
 * The control flags can be a combination (bitwise OR) of the following:
 *
 *   - ZRAN_INFLATE_INIT_Z_STREAM:         Initialise the z_stream struct
 *                                         before inflation.
 *
 *   - ZRAN_INFLATE_FREE_Z_STREAM:         Clean up the z_stream struct
 *                                         after inflation.
 *
 *   - ZRAN_INFLATE_INIT_READBUF:          Allocate a read buffer before
 *                                         inflation.
 *
 *   - ZRAN_INFLATE_FREE_READBUF:          Free the read buffer after
 *                                         inflation.
 *
 *   - ZRAN_INFLATE_USE_OFFSET:            If set, use the provided offset
 *                                         parameter; otherwise, use the
 *                                         offsets stored in the index
 *                                         struct.
 *
 *   - ZRAN_INFLATE_CLEAR_READBUF_OFFSETS: If set, clear the read buffer
 *                                         offset/length stored in the index
 *                                         struct, otherwise assume that they
 *                                         are valid.
 *
 *   - ZRAN_INFLATE_STOP_AT_BLOCK:         If set, stop inflation when a
 *                                         deflate block boundary is reached
 *                                         (the Z_BLOCK flag is passed to the
 *                                         zlib inflate function). Otherwise,
 *                                         inflation will continue until one
 *                                         of the conditions in step 7, above,
 *                                         are met.
 *
 * This function returns one of the following codes. Furthermore, if an error
 * did not occur (i.e. anything but ZRAN_INFLATE_ERROR or
 * ZRAN_INFLATE_NOT_COVERED was returned), the total_consumed and total_output
 * parameters are respectively updated to contain the total number of
 * compressed bytes that were read from the file, and total number of
 * decompressed bytes that were copied to the data buffer.
 *
 *   - ZRAN_INFLATE_OK:             Inflation was successful and the requested
 *                                  number of bytes were copied to the provided
 *                                  data buffer.
 *
 *   - ZRAN_INFLATE_NOT_COVERED:    The requested compressed data offset is not
 *                                  covered by the index.
 *
 *   - ZRAN_INFLATE_OUTPUT_FULL:    The provided data buffer has been filled.
 *
 *   - ZRAN_INFLATE_BLOCK_BOUNDARY: A deflate block boundary was encountered.
 *                                  This will only be returned if the
 *                                  ZRAN_INFLATE_STOP_AT_BLOCK flag is active.
 *
 *   - ZRAN_INFLATE_EOF:            The end of file has been reached.
 *
 *   - ZRAN_INFLATE_CRC_ERROR:      The CRC or uncompressed data size in the
 *                                  GZIP footer does not match the CRC/size
 *                                  that was calculated.
 *
 *   - ZRAN_INFLATE_ERROR:          A critical error has occurred.
 */
static int _zran_inflate(
    zran_index_t *index,            /* Pointer to the index. */

    z_stream     *strm,             /* Pointer to a z_stream struct. */

    uint64_t      offset,           /* Compressed data offset to start
                                       inflation from. */

    uint16_t      flags,            /* Control flags. */

    uint32_t     *total_consumed,   /* Pointer which is updated to contain the
                                       total number of bytes that were read
                                       from the input file. */

    uint32_t     *total_output,     /* Pointer which is updated to contain the
                                       total number of bytes that were
                                       inflated, and stored in data. */

    uint32_t      len,              /* Maximum number of bytes to inflate. May
                                       be 0. */

    uint8_t      *data,             /* Place to store the inflated bytes. */

    int           add_stream_points /* Add index points at the beginning of
                                       every gzip stream, including the first
                                       one at the beginning of the input file
                                     */
);


/* Initialise a zran_index_t struct for use with the given GZIP file. */
int zran_init(zran_index_t *index,
              FILE         *fd,
              PyObject     *f,
              uint32_t      spacing,
              uint32_t      window_size,
              uint32_t      readbuf_size,
              uint16_t      flags)
{

    zran_point_t *point_list = NULL;
    int64_t       compressed_size;

    zran_log("zran_init(%u, %u, %u, %u)\n",
             spacing, window_size, readbuf_size, flags);

    if (spacing      == 0) spacing      = 1048576;
    if (window_size  == 0) window_size  = 32768;
    if (readbuf_size == 0) readbuf_size = 16384;

    /*
     * The zlib manual specifies that a window
     * size of 32KB is 'always enough' to
     * initialise inflation/deflation with a
     * set dictionary. Less than that is not
     * guaranteed to be enough.
     */
    if (window_size < 32768)
        goto fail;

    /*
     * Small read-buffers make code complicated.
     * The absolute minimum we need is enough to
     * store a GZIP footer, null padding bytes at
     * the end of a stream, and the subsequent
     * GZIP header. There are no bounds on the
     * number of padding bytes, or the size of a
     * GZIP header, so this constraint is
     * arbitrary (but should be good enough).
     */
    if (readbuf_size < 128)
        goto fail;

    /*
     * window_size bytes of uncompressed data are
     * stored with each seek point in the index.
     * So it's a bit silly to have the distance
     * between consecutive points less than the
     * window size.
     */
    if (spacing <= window_size)
        goto fail;

    /* The file must be opened in read-only mode */
    if (!is_readonly(fd, f))
        goto fail;

    /*
     * Calculate the size of the compressed file
     */
    if (seekable_(fd, f)) {
        if (fseek_(fd, f, 0, SEEK_END) != 0)
            goto fail;

        compressed_size = ftell_(fd, f);

        if (compressed_size < 0)
            goto fail;

        if (fseek_(fd, f, 0, SEEK_SET) != 0)
            goto fail;
    } else {
        /*
         * File is not seekable, so don't calculate
         * compressed_size.  It will be updated in
         * _zran_read_data_from_file when the EOF
         * is reached
         */
        compressed_size = 0;
    }

    /*
     * Allocate some initial space
     * for the index point list
     */
    point_list = calloc(1, sizeof(zran_point_t) * 8);
    if (point_list == NULL) {
        goto fail;
    }

    /* initialise the index struct */
    index->fd                   = fd;
    index->f                    = f;
    index->flags                = flags;
    index->compressed_size      = compressed_size;
    index->uncompressed_size    = 0;
    index->spacing              = spacing;
    index->window_size          = window_size;
    index->log_window_size      = (int)round(log10(window_size) / log10(2));
    index->readbuf_size         = readbuf_size;
    index->readbuf_offset       = 0;
    index->readbuf_end          = 0;
    index->readbuf              = NULL;
    index->npoints              = 0;
    index->size                 = 8;
    index->uncmp_seek_offset    = 0;
    index->inflate_cmp_offset   = 0;
    index->inflate_uncmp_offset = 0;
    index->validating           = 0;
    index->last_stream_ended    = 0;
    index->stream_size          = 0;
    index->stream_crc32         = 0;
    index->list                 = point_list;

    return 0;

fail:
    free(point_list);
    return -1;
}


/* Returns the compressed or uncompressed index limit. */
uint64_t _zran_index_limit(zran_index_t *index, uint8_t compressed) {

    if (index->npoints == 0)
        return 0;

    if (compressed) return index->list[index->npoints - 1].cmp_offset;
    else            return index->list[index->npoints - 1].uncmp_offset;
}


/* Expands the memory used to store the index points. */
int _zran_expand_point_list(zran_index_t *index) {
    zran_point_t *new_list;
    uint32_t new_size = index->size * 2;

    zran_log("_zran_expand_point_list(%i -> %i)\n", index->size, new_size);

    new_list = realloc(index->list, sizeof(zran_point_t) * new_size);

    if (new_list == NULL) {
        /* old list is still valid */
        return -1;
    }

    index->list = new_list;
    index->size = new_size;

    return 0;
}


/* Frees any unused memory allocated for index storage. */
int _zran_free_unused(zran_index_t *index) {

    zran_point_t *new_list;
    size_t        new_size;

    zran_log("_zran_free_unused\n");

    if (index->npoints < 8) new_size = 8;
    else                    new_size = index->npoints;

    new_list = realloc(index->list, sizeof(zran_point_t) * new_size);

    if (new_list == NULL) {
        return -1;
    }

    index->list = new_list;
    index->size = new_size;

    return 0;
}


/* Deallocate memory used by a zran_index_t struct. */
void zran_free(zran_index_t *index) {

    uint32_t      i;
    zran_point_t *pt;

    zran_log("zran_free\n");

    for (i = 0; i < index->npoints; i++) {
        pt = &(index->list[i]);

        /*
         * points at compression stream boundaries
         * have no data associated with them
         */
        if (pt->data != NULL) {
            free(pt->data);
        }
    }

    free(index->list);

    index->fd                = NULL;
    index->f                 = NULL;
    index->spacing           = 0;
    index->window_size       = 0;
    index->readbuf_size      = 0;
    index->npoints           = 0;
    index->size              = 0;
    index->list              = NULL;
    index->uncmp_seek_offset = 0;
}


/* Discard all points in the index after the specified compressed offset. */
int _zran_invalidate_index(zran_index_t *index, uint64_t from)
{
    uint64_t      i;
    zran_point_t *p;

    if (index->npoints == 0)
        return 0;

    for (i = 0; i < index->npoints; i++) {

        p = &(index->list[i]);

        if (p->cmp_offset >= from)
            break;
    }

    /*
     * The index doesn't cover
     * the requested offest
     */
    if (i == index->npoints)
        return 0;

    if (i <= 1) index->npoints = 0;
    else        index->npoints = i - 1;

    return _zran_free_unused(index);
}


/* (Re-)Builds the full index. */
int zran_build_index(zran_index_t *index, uint64_t from, uint64_t until)
{

    if (_zran_invalidate_index(index, from) != 0)
        return ZRAN_BUILD_INDEX_FAIL;

    return _zran_expand_index(index, until);
}


/* Searches for and returns the index at the specified offset. */
int _zran_get_point_at(
    zran_index_t  *index,
    uint64_t       offset,
    uint8_t        compressed,
    zran_point_t **point)
{
    uint64_t      cmp_max;
    uint64_t      uncmp_max;
    zran_point_t *last;
    zran_point_t *prev;
    zran_point_t *curr;
    uint8_t       bit;
    uint32_t      i;

    *point = NULL;

    /*
     * Bad input - past the end of the compressed or
     * uncompressed streams (if their sizes are known).
     */
    if (compressed                 &&
        index->compressed_size > 0 &&
        offset >= index->compressed_size)
        goto eof;

    if (!compressed                  &&
        index->uncompressed_size > 0 &&
        offset >= index->uncompressed_size)
        goto eof;

    if (index->npoints == 0)
        goto not_covered;

    zran_log("_zran_get_point_at(%llu, c=%u)\n", offset, compressed);

    /*
     * Figure out how much of the compressed
     * and uncompressed data the index currently
     * covers -  the offsets of the last point
     * in the index.
     */
    last      = &(index->list[index->npoints - 1]);
    uncmp_max = last->uncmp_offset;
    cmp_max   = last->cmp_offset;

    if ( compressed && offset > cmp_max)   goto not_covered;
    if (!compressed && offset > uncmp_max) goto not_covered;

    /*
     * We should have an index point
     * which corresponds to this offset,
     * so let's search for it.
     */
    prev = index->list;
    for (i = 1; i < index->npoints; i++) {

        curr = &(index->list[i]);

        if (compressed) {

            /*
             * Adjust the offset for non
             * byte-aligned seek points.
             */
            if (curr->bits > 0) bit = 1;
            else                bit = 0;

            if (curr->cmp_offset > offset + bit)
                break;

        }
        else {
            if (curr->uncmp_offset > offset)
                break;
        }

        prev = curr;
    }

    *point = prev;
    return ZRAN_GET_POINT_OK;

not_covered:
    *point = NULL;
    return ZRAN_GET_POINT_NOT_COVERED;

eof:
    *point = NULL;
    return ZRAN_GET_POINT_EOF;
}


/*
 * Get the index point corresponding to the given offset, expanding
 * the index as needed if ZRAN_AUTO_BUILD is active.
 */
int _zran_get_point_with_expand(zran_index_t  *index,
                                uint64_t       offset,
                                uint8_t        compressed,
                                zran_point_t **point)
{

    int      result;
    uint64_t expand;
    uint64_t limit;

    zran_log("_zran_get_point_with_expand(%llu, %u, autobuild=%u)\n",
             offset,
             compressed,
             index->flags & ZRAN_AUTO_BUILD);

    /*
     * See if there is an index point that
     * covers the specified offset. If there's
     * not, we're going to expand the index
     * until there is.
     */
    result = _zran_get_point_at(index, offset, compressed, point);

    /*
     * Don't expand the index if
     * auto_build is not active
     */
    if ((index->flags & ZRAN_AUTO_BUILD) == 0) {
        return result;
    }

    while (result == ZRAN_GET_POINT_NOT_COVERED) {

        /*
         * If result == ZRAN_GET_POINT_NOT_COVERED,
         * get_point says that an index point for
         * this offset doesn't yet exist. So let's
         * expand the index.
         *
         * Guess how far we need to expand the index,
         * and expand it by that much.
         */
        if (compressed == 0) expand = _zran_estimate_offset(index, offset, 0);
        else                 expand = offset;

        /*
         * If _zran_estimate_offset was unable to
         * estimate a sensible compressed offset
         * (i.e. smaller or at the current index
         * extent), we force it past the limit,
         * so that the expand_index function will
         * create at least one point.
         */
        limit = _zran_index_limit(index, 1);
        if (expand <= limit)
            expand = limit + 10;

        zran_log("Estimated mapping from uncompressed offset "
                 "%lu into compressed data: %lu\n", offset, expand);

        /*
         * Expand the index
         */
        result = _zran_expand_index(index, expand);
        if      (result == ZRAN_EXPAND_INDEX_CRC_ERROR) { goto crcerror; }
        else if (result != 0)                           { goto fail; }

        /*
         * Index has been expanded, so
         * there should now be a point
         * which covers the requested
         * offset.
         */
        result = _zran_get_point_at(index, offset, compressed, point);

        /*
         * If we've made it to EOF, return
         * a ref to the eof point.
         */
        if (result == ZRAN_GET_POINT_EOF) {
            *point = &index->list[index->npoints - 1];

            if (offset < index->uncompressed_size) {
                result = ZRAN_GET_POINT_OK;
            }
        }
    }

    return result;

crcerror:
    return ZRAN_GET_POINT_CRC_ERROR;
fail:
    return ZRAN_GET_POINT_FAIL;
}


/*
 * Given an offset in one stream, estimates the corresponding offset into the
 * other stream.
 */
uint64_t _zran_estimate_offset(
    zran_index_t *index,
    uint64_t      offset,
    uint8_t       compressed)
{

    zran_point_t *last;
    uint64_t      estimate;

    /*
     * The first index in the list maps
     * indices 0 and 0, which won't help
     * us here. So we need at least two
     * index points.
     */
    if (index->npoints <= 1) last = NULL;
    else                     last = &(index->list[index->npoints - 1]);

    /*
     * We have no reference. At least two
     * index points need to have been created.
     * The assumed correspondences between
     * the compressed streams are arbitrary.
     */
    if (last == NULL) {
        if (compressed) estimate = offset * 2.0;
        else            estimate = offset * 0.8;
    }

    /*
     * I'm just assuming a roughly linear correspondence
     * between the compressed/uncompressed data streams.
     */
    else if (compressed) {
        estimate = round(offset *
                         ((float)last->uncmp_offset / last->cmp_offset));
    }
    else {
        estimate = round(offset *
                         ((float)last->cmp_offset / last->uncmp_offset));
    }


    zran_log("_zran_estimate_offset(%llu, %u) = %llu\n",
             offset, compressed, estimate);

    return estimate;
}



/* Add a new point to the index. */
int _zran_add_point(zran_index_t  *index,
                    uint8_t        bits,
                    uint64_t       cmp_offset,
                    uint64_t       uncmp_offset,
                    uint32_t       data_offset,
                    uint32_t       data_size,
                    uint8_t       *data) {

    uint8_t      *point_data = NULL;
    zran_point_t *next       = NULL;

    #ifdef ZRAN_VERBOSE
    zran_log("_zran_add_point(%i, c=%lld + %i, u=%lld, data=%u / %u)\n",
             index->npoints,
             cmp_offset,
             bits > 0,
             uncmp_offset,
             data_offset,
             data_size);

    if (data != NULL)
        zran_log("Window data: [%02x %02x %02x %02x ...]\n",
                 data[(data_offset - index->window_size + 0) % data_size],
                 data[(data_offset - index->window_size + 1) % data_size],
                 data[(data_offset - index->window_size + 2) % data_size],
                 data[(data_offset - index->window_size + 3) % data_size]);
    #endif

    /* if list is full, make it bigger */
    if (index->npoints == index->size) {
        if (_zran_expand_point_list(index) != 0) {
            goto fail;
        }
    }

    /*
     * Allocate memory to store the
     * uncompressed data (the "window")
     * associated with this point. Index
     * points corresponding to the beginning
     * of a gzip stream (including at start
     * of file) do not have any window data
     * associated with them.
     */
    if (data == NULL) {
        point_data = NULL;
    }
    else {
        point_data = calloc(1, index->window_size);
        if (point_data == NULL)
            goto fail;
    }

    next               = &(index->list[index->npoints]);
    next->bits         = bits;
    next->cmp_offset   = cmp_offset;
    next->uncmp_offset = uncmp_offset;
    next->data         = point_data;

    /*
     * The uncompressed data may not start at
     * the beginning of the data pointer, but
     * rather from an arbitrary point. So we
     * copy the beginning of the window from
     * the end of data, and the end of the
     * window from the beginning of data. Does
     * that make sense?
     */
    if (data != NULL) {
        if (data_offset >= index->window_size) {

            memcpy(point_data,
                   data + (data_offset - index->window_size),
                   index->window_size);

            zran_log("Copy %u bytes from %u to %u\n",
                     index->window_size,
                     data_offset - index->window_size,
                     data_offset);
        }
        else {
            memcpy(point_data,
                   data + (data_size - (index->window_size - data_offset)),
                   (index->window_size - data_offset));

            memcpy(point_data + (index->window_size - data_offset),
                   data,
                   data_offset);

            zran_log("Copy %u bytes from %u to %u, %u bytes from %u to %u\n",
                     (index->window_size - data_offset),
                     (data_size - (index->window_size - data_offset)),
                     data_size,
                     data_offset,
                     0,
                     data_offset);
        }
    }

    index->npoints++;

    return 0;

fail:
    free(point_data);

    return -1;
}


/* Initialise the given z_stream struct for decompression/inflation. */
int _zran_init_zlib_inflate(zran_index_t *index,
                            z_stream     *strm,
                            zran_point_t *point) {

    int           ret;
    int           window;
    int64_t       seek_loc;
    unsigned long bytes_read;

    bytes_read   = strm->avail_in;
    window       = index->log_window_size;
    strm->zalloc = Z_NULL;
    strm->zfree  = Z_NULL;
    strm->opaque = Z_NULL;

    /*
     * If we're starting from the the current location in
     * the compressed data, we assume that there is a gzip
     * header present. Initialise inflation, then read
     * past the header.

     * Below, we will re-initialise for raw inflation.
     */
    if (point == NULL) {

        zran_log("_zran_init_zlib_inflate from current "
                 "seek location (expecting GZIP header)\n");
        if (inflateInit2(strm, window + 32) != Z_OK) { goto fail; }
        if (inflate(strm, Z_BLOCK)          != Z_OK) { goto fail_free_strm; }
        if (inflateEnd(strm)                != Z_OK) { goto fail; }
    }

    /*
     * If starting from an index point, seek to the
     * required location in the compressed data stream.
     *
     * The compressed offset for index points correspond
     * to the first full byte of compressed data. So if
     * the index point is not byte-aligned (bits > 0), we
     * need to seek to the previous byte, and tell zlib
     * about it (via the inflatePrime call below).
     */
    else {
        seek_loc = point->cmp_offset - (point->bits > 0);

        zran_log("_zran_init_zlib_inflate from index point (%li, %li, %li)\n",
                 seek_loc,
                 point->cmp_offset,
                 point->uncmp_offset);

        if (fseek_(index->fd, index->f, seek_loc, SEEK_SET) != 0) {
            goto fail;
        }
    }

    /*
     * Now initialise for raw inflation. This tells zlib
     * not to expect any GZIP headers, and not to read
     * the GZIP footer (as we do our own CRC validation
     * in _zran_inflate).
     */

    if (inflateInit2(strm, -window) != Z_OK) {
        goto fail;
    }

    /*
     * If starting from an index point, initialise
     * the inflation dictionary from the uncompressed
     * data associated with the index point.
     */
    if (point != NULL && point->data != NULL) {

        /*
         * The starting index point is not byte-aligned,
         * so we'll insert the initial bits into the
         * inflate stream using inflatePrime (above,
         * we seeked one byte back to accommodate this).
         */
        if (point->bits > 0) {

            ret = getc_(index->fd, index->f);

            if (ret == -1 && ferror_(index->fd, index->f)) {
                goto fail_free_strm;
            }

            if (inflatePrime(strm,
                             point->bits, ret >> (8 - point->bits)) != Z_OK)
                goto fail_free_strm;
        }

        /*
         * Initialise the inflate stream
         * with the index point data.
         */
        if (inflateSetDictionary(strm,
                                 point->data,
                                 index->window_size) != Z_OK)
            goto fail_free_strm;
    }

    /*
     * Reset CRC/size validation counters when
     * we start reading a new gzip stream
     */
    index->validating   = (point == NULL);
    index->stream_size  = 0;
    index->stream_crc32 = 0;

    zran_log("_zran_zlib_init_inflate: initialised, read %i bytes\n",
             bytes_read - strm->avail_in);

    /*
     * Return the number of bytes of compressed
     * data, if any that were read over
     */
    return bytes_read - strm->avail_in;

/*
 * Something has gone wrong, but
 * inflateInit2 has been called,
 * so we need to call inflateEnd.
 * Falls through to the fail:
 * clause.
 */
fail_free_strm:
    inflateEnd(strm);
/* Something has gone wrong */
fail:
    return -1;
}


/*
 * Read data from the GZIP file, and copy it into the read buffer for
 * decompression.
 */
static int _zran_read_data_from_file(zran_index_t *index,
                                     z_stream     *stream,
                                     uint64_t      cmp_offset,
                                     uint64_t      uncmp_offset,
                                     uint32_t      need_atleast) {

    size_t f_ret;

    if (stream->avail_in >= need_atleast) {
        return 0;
    }

    /*
     * If there are any unprocessed bytes
     * left over, put them at the beginning
     * of the read buffer.
     *
     * TODO: In times gone by, we would only
     * attempt to read data (and therefore
     * rotate memory here) when the read
     * buffer was empty. But now, to keep
     * the code in _zran_inflate clean-ish,
     * we do this repeatedly, even when we
     * are at EOF, to ensure that there is
     * enough data to validate one stream,
     * and find the next. We could improve
     * things here, by only rotating memory
     * here if needed.
     */
    if (stream->avail_in > 0) {
        memmove(index->readbuf, stream->next_in, stream->avail_in);
    }

    zran_log("Reading from file %llu "
             "[into readbuf offset %u]\n",
             cmp_offset + stream->avail_in,
             stream->avail_in);

    /*
     * Read a block of compressed data
     * (offsetting past any left over
     * bytes that we may have copied to
     * the beginning of the read buffer
     * above).
     */
    f_ret = fread_(index->readbuf + stream->avail_in,
                   1,
                   index->readbuf_size - stream->avail_in,
                   index->fd,
                   index->f);

    if (ferror_(index->fd, index->f)) {
        goto fail;
    }

    /*
     * No bytes left to read, and there are
     * only 8 bytes left to process (size of
     * gzip footer) - we've reached EOF.
     */
    if (f_ret == 0 && stream->avail_in <= 8) {
        if (feof_(index->fd, index->f, f_ret)) {

            zran_log("End of file, stopping inflation\n");

            /*
             * Reset next_in pointer to beginning of
             * read buffer, as we rotated it above,
             * and the area that next_in was pointing
             * to may have been overwritten by memmove.
             */
            stream->next_in = index->readbuf;

            /*
             * we have uncompressed everything,
             * so we now know its size.
             */
            if (index->uncompressed_size == 0) {
                zran_log("Updating uncompressed data "
                         "size: %llu\n", uncmp_offset);
                index->uncompressed_size = uncmp_offset;
            }
            if (index->compressed_size == 0) {
                zran_log("Updating compressed data "
                         "size: %llu\n", cmp_offset);
                index->compressed_size = cmp_offset + 8;
            }
            goto eof;
        }
        /*
         * Or something went wrong (this
         * should never happen if ferror
         * does the right thing).
         */
        else {
            goto fail;
        }
    }

    zran_log("Read %lu bytes from file [c=%llu, u=%llu] "
             "[%02x %02x %02x %02x ...]\n",
             f_ret, cmp_offset, uncmp_offset,
             index->readbuf[stream->avail_in],
             index->readbuf[stream->avail_in + 1],
             index->readbuf[stream->avail_in + 2],
             index->readbuf[stream->avail_in + 3]);

    /*
     * Tell zlib about the block
     * of compressed data that we
     * just read in.
     */
    index->readbuf_end = f_ret + stream->avail_in;
    stream->avail_in  += f_ret;
    stream->next_in    = index->readbuf;

    return 0;
eof:
    return ZRAN_READ_DATA_EOF;
fail:
    return ZRAN_READ_DATA_ERROR;
}


/*
 * Identify the location of the next compressed stream (if the file
 * contains concatenated streams).
 */
int _zran_find_next_stream(zran_index_t *index,
                           z_stream     *stream,
                           int          *offset) {

    int ret;
    int found;

    /*
     * Search for the beginning of
     * the next stream. GZIP files
     * start with 0x1f8b.
     */
    found = 0;

    zran_log("Searching for a new stream [%u]\n", stream->avail_in);

    while (stream->avail_in > 0) {

        if (stream->avail_in   >= 2    &&
            stream->next_in[0] == 0x1f &&
            stream->next_in[1] == 0x8b) {
            found = 1;
            break;
        }

        *offset          += 1;
        stream->next_in  += 1;
        stream->avail_in -= 1;
    }

    /*
     * No header found for
     * the next stream.
     */
    if (found == 0) {
        zran_log("Could not find another stream [%u]\n", stream->avail_in);
        goto not_found;
    }

    zran_log("New stream found, re-initialising inflation\n");

    /*
     * Re-configure for inflation
     * from the new stream.
     */
    if (inflateEnd(stream) != Z_OK) {
        goto fail;
    }

    ret = _zran_init_zlib_inflate(index, stream, NULL);

    if (ret < 0) {
        goto fail;
    }

    *offset += ret;

    return 0;

fail:
    return ZRAN_FIND_STREAM_ERROR;

not_found:
    return ZRAN_FIND_STREAM_NOT_FOUND;
}


/* Validate the CRC32 and size of a GZIP stream. */
static int _zran_validate_stream(zran_index_t *index,
                                 z_stream     *stream,
                                 int          *offset) {

    uint32_t crc;
    uint32_t size;

    /* CRC validation is disabled */
    if (index->flags & ZRAN_SKIP_CRC_CHECK) {
        return 0;
    }

    /*
     * A gzip stream should end with an 8 byte footer,
     * which contains the CRC32 of the uncompressed
     * data, and the uncompressed size modulo 2^32.
     */
    if (stream->avail_in < 8) {
        return ZRAN_VALIDATE_STREAM_ERROR;
    }

    crc  = ((stream->next_in[0] << 0)  +
            (stream->next_in[1] << 8)  +
            (stream->next_in[2] << 16) +
            (stream->next_in[3] << 24));
    size = ((stream->next_in[4] << 0)  +
            (stream->next_in[5] << 8)  +
            (stream->next_in[6] << 16) +
            (stream->next_in[7] << 24));

    zran_log("Validating CRC32 and size [%8x == %8x, %u == %u]\n",
             crc, index->stream_crc32, size, index->stream_size);

    stream->avail_in -= 8;
    stream->next_in  += 8;
    *offset          += 8;

    if (index->stream_crc32 != crc || index->stream_size != size) {
        return ZRAN_VALIDATE_STREAM_INVALID;
    }

    return 0;
}


/* The workhorse. Inflate/decompress data from the file. */
static int _zran_inflate(zran_index_t *index,
                         z_stream     *strm,
                         uint64_t      offset,
                         uint16_t      flags,
                         uint32_t     *total_consumed,
                         uint32_t     *total_output,
                         uint32_t      len,
                         uint8_t      *data,
                         int           add_stream_points) {

    /*
     * z_ret is for zlib/zran functions.
     * off is for storing offset of new
     * stream (from _zran_validate_stream
     * and _zran_find_next_stream).
     * return_val/error_return_val is
     * the return value for this function.
     */
    int z_ret;
    int off;
    int return_val       = ZRAN_INFLATE_OK;
    int error_return_val = ZRAN_INFLATE_ERROR;

    /*
     * Offsets into the compressed
     * and uncompressed data streams,
     * and total number of bytes
     * decompressed and output.
     */
    uint64_t cmp_offset;
    uint64_t uncmp_offset;
    uint32_t _total_consumed = 0;
    uint32_t _total_output   = 0;

    /*
     * Number of bytes input/decompressed
     * during a single call to zlib:inflate.
     */
    uint32_t bytes_consumed = 0;
    uint32_t bytes_output   = 0;

    /*
     * Index point to start from
     * (if ZRAN_INFLATE_USE_OFFSET
     * is active).
     */
    zran_point_t *start = NULL;

    /*
     * Set all zstream_t fields to 0
     * if we are initialising.
     */
    if (inflate_init_stream(flags)) {
        memset(strm, 0, sizeof(z_stream));
    }

    /*
     * If ZRAN_INFLATE_INIT_READBUF is not set,
     * make sure that a read buffer exists.
     *
     * If the opposite is true, the read buffer
     * from a prior call has not been cleaned up.
     */
    if ((!inflate_init_readbuf(flags) && index->readbuf == NULL) ||
        ( inflate_init_readbuf(flags) && index->readbuf != NULL)) {
        goto fail;
    }

    /*
     * It begins...
     */
    zran_log("_zran_inflate(%llu, block=%u, use_offset=%u, init_stream=%u,\n"
             "              free_stream=%u, init_readbuf=%u, free_readbuf=%u,\n"
             "              clear_offsets=%u, nbytes=%u)\n",
             offset,
             inflate_stop_at_block(        flags),
             inflate_use_offset(           flags),
             inflate_init_stream(          flags),
             inflate_free_stream(          flags),
             inflate_init_readbuf(         flags),
             inflate_free_readbuf(         flags),
             inflate_clear_readbuf_offsets(flags),
             len);

    /*
     * The compressed/uncompressed offsets are initialised in
     * one of three ways. If ZRAN_INFLATE_USE_OFFSET is active,
     * they are either:
     *
     *    - Both set to 0
     *
     *    - Initialised according to an existing index
     *      point that precedes the requested offset.
     *
     * Otherwise, they are initialised from index->inflate_cmp_offset
     * and index->inflate_uncmp_offset, which are assumed to have been
     * set in a prior call to _zran_inflate.
     */
    if (inflate_use_offset(flags)) {

        cmp_offset   = 0;
        uncmp_offset = 0;

        /*
         * If a non-zero offset has been specified,
         * search the index to see if we can start
         * inflating from a known location.
         */
        if (offset > 0) {

            /*
             * In order to successfully decompress
             * data from the current uncompressed seek
             * location, we need to start decompressing
             * from the index point which preceeds it.
             */
            z_ret = _zran_get_point_at(index, offset, 1, &start);

            if (z_ret == ZRAN_GET_POINT_NOT_COVERED)
                return ZRAN_INFLATE_NOT_COVERED;

            if (z_ret == ZRAN_GET_POINT_EOF)
                return ZRAN_INFLATE_EOF;
        }

        /*
         * Start inflating from the index point
         * corresponding to the offset (or keep
         * the offsets at 0 if no point was found).
         */
        if (start != NULL) {

            cmp_offset   = start->cmp_offset;
            uncmp_offset = start->uncmp_offset;
        }
    }

    /*
     * If ZRAN_INFLATE_USE_OFFSET is not active,
     * we initialise from offsets which were
     * stored on the last call to _zran_inflate.
     */
    else {
        cmp_offset   = index->inflate_cmp_offset;
        uncmp_offset = index->inflate_uncmp_offset;
    }

    zran_log("initialising to inflate from "
             "cmp_offset=%llu, uncmp_offset=%llu\n",
             cmp_offset,
             uncmp_offset);

    /*
     * If ZRAN_INFLATE_INIT_READBUF,
     * allocate memory for reading
     * compressed data from the file.
     * The buffer is attached to the
     * zran_index_t->readbuf pointer.
     */
    if (inflate_init_readbuf(flags)) {
        index->readbuf = calloc(1, index->readbuf_size);
        if (index->readbuf == NULL)
            goto fail;
    }

    /*
     * If ZRAN_INFLATE_CLEAR_READBUF_OFFSETS,
     * we clear any stored information about
     * the read buffer, and start reading
     * from/writing to it from the beginning.
     */
    if (inflate_clear_readbuf_offsets(flags)) {
        index->readbuf_offset = 0;
        index->readbuf_end    = 0;
    }

    /*
     * Otherwise, assume that there is already
     * some input (compressed) data in the
     * readbuf, and that index->readbuf_offset
     * and index->readbuf_end were set on a
     * prior call.
     *
     *    - readbuf_offset tells us where in
     *      readbuf the data starts
     *
     *    - readbuf_end tells us where it ends.
     */
    else {
        strm->next_in  = index->readbuf     + index->readbuf_offset;
        strm->avail_in = index->readbuf_end - index->readbuf_offset;
    }

    /*
     * Tell zlib where to store
     * the uncompressed data.
     */
    strm->avail_out = len;
    strm->next_out  = data;

    /*
     * If ZRAN_INFLATE_INIT_Z_STREAM is active,
     * initialise the zlib struct for inflation.
     *
     * If ZRAN_INFLATE_INIT_Z_STREAM is not
     * active, we assume that the inflation is
     * already initialised.
     */
    if (inflate_init_stream(flags)) {

        /*
         * No index point - we need to start reading
         * from the beginning of the input file. If
         * starting from an index point, the
         * _zran_init_zlib_inflate function will seek
         * to the correct location in the file for us.
         */
        if (start == NULL) {
            /*
             * If file is not seekable, assume that
             * it is positioned at the beginning of
             * the stream.
             */
            if (seekable_(index->fd, index->f)) {
                if (fseek_(index->fd, index->f, 0, SEEK_SET) != 0) {
                    goto fail;
                }
            }

            /*
             * In this situation, _zran_init_zlib_inflate
             * is going to expect a GZIP header, so make
             * sure we have some data for it to look at.
             */
            if (_zran_read_data_from_file(index,
                                          strm,
                                          cmp_offset,
                                          uncmp_offset,
                                          index->readbuf_size) != 0) {
                goto fail;
            }
        }

        /*
         * If init_zlib_inflate skips over any input data
         * (e.g. gzip header), it returns the number of
         * bytes that were read.
         */
        z_ret = _zran_init_zlib_inflate(index, strm, start);
        if (z_ret < 0) {
            goto fail;
        }
        cmp_offset      += z_ret;
        _total_consumed += z_ret;

        if (start == NULL && add_stream_points) {
            if (_zran_add_point(index, 0, cmp_offset, 0, 0, 0, NULL) != 0) {
                goto fail;
            }
        }
    }

    /*
     * Keep going until we run out of space.
     */
    while (strm->avail_out > 0) {

        /*
         * Make sure the input buffer contains
         * some data to be decompressed.
         */
        z_ret = _zran_read_data_from_file(index,
                                          strm,
                                          cmp_offset,
                                          uncmp_offset,
                                          1);

        /* EOF - there is no more data left to read */
        if (z_ret == ZRAN_READ_DATA_EOF) {
            return_val = ZRAN_INFLATE_EOF;
            break;
        }
        else if (z_ret != 0) {
            goto fail;
        }

        /*
         * Decompress data until there's no data
         * left, or we've read enough bytes
         */
        z_ret = Z_OK;
        while (strm->avail_in > 0) {

            /*
             * Initialise counters to calculate
             * how many bytes are input/output
             * during this call to inflate.
             */
            bytes_consumed = strm->avail_in;
            bytes_output   = strm->avail_out;

            zran_log("Before inflate - avail_in=%u, avail_out=%u, "
                     "cmp_offset=%lu, uncmp_offset=%lu\n",
                     strm->avail_in, strm->avail_out,
                     cmp_offset, uncmp_offset);

            /*
             * Inflate the block - the decompressed
             * data is output straight to the provided
             * data buffer.
             *
             * If ZRAN_INFLATE_STOP_AT_BLOCK is active,
             * Z_BLOCK tells inflate to stop inflating
             * at a compression block boundary. Otherwise,
             * inflate stops when it comes to the end of a
             * stream, or it runs out of input or output.
             */
            if (inflate_stop_at_block(flags)) {
                z_ret = inflate(strm, Z_BLOCK);
            }
            else {
                z_ret = inflate(strm, Z_NO_FLUSH);
            }

            /*
             * Adjust our offsets according to what
             * was actually consumed/decompressed.
             */
            bytes_consumed   = bytes_consumed - strm->avail_in;
            bytes_output     = bytes_output   - strm->avail_out;
            cmp_offset      += bytes_consumed;
            _total_consumed += bytes_consumed;
            uncmp_offset    += bytes_output;
            _total_output   += bytes_output;

            zran_log("After inflate - avail_in=%u, avail_out=%u, "
                     "cmp_offset=%lu, uncmp_offset=%lu\n",
                     strm->avail_in, strm->avail_out,
                     cmp_offset, uncmp_offset);


            /*
             * Now we need to figure out what just happened.
             *
             * Z_BUF_ERROR indicates that the output buffer
             * is full; we clobber it though, as it makes the
             * code below a bit easier (and anyway, we can
             * tell if the output buffer is full by checking
             * strm->avail_out).
             */
            if (z_ret == Z_BUF_ERROR) {
                z_ret = Z_OK;
            }

            /*
             * If z_ret is not Z_STREAM_END or
             * Z_OK, something has gone wrong.
             *
             * If the file comprises a sequence of
             * concatenated gzip streams, we will
             * encounter Z_STREAM_END before the end
             * of the file (where one stream ends and
             * the other begins).
             *
             * If at a new stream, we re-initialise
             * inflation on the next loop iteration.
             */
            if (z_ret != Z_OK && z_ret != Z_STREAM_END) {
                zran_log("zlib inflate failed (code: %i, msg: %s)\n",
                         z_ret, strm->msg);
                goto fail;
            }

            /*
             * If we have not yet validated the current
             * GZIP stream, update its size and crc so
             * we can validate them against the size
             * recorded in the stream footer when we
             * get to it.
             */
            if ((uncmp_offset > index->last_stream_ended) &&
                index->validating                         &&
                !(index->flags & ZRAN_SKIP_CRC_CHECK)) {
                index->stream_size +=       bytes_output;
                index->stream_crc32 = crc32(index->stream_crc32,
                                            strm->next_out - bytes_output,
                                            bytes_output);
            }

            /*
             * End of a block? If INFLATE_STOP_AT_BLOCK
             * is active, we want to stop at a compression
             * block boundary.
             *
             * If we used Z_BLOCK above, and inflate
             * encountered a block boundary, it indicates
             * this in the the strm->data_type field.
             */
            if (inflate_stop_at_block(flags) &&
                ((strm->data_type & 128) && !(strm->data_type & 64))) {

                zran_log("At block or stream boundary, "
                         "stopping inflation\n");

                return_val = ZRAN_INFLATE_BLOCK_BOUNDARY;
                break;
            }

            /*
             * We've found the end of file, or end of one
             * gzip stream. Validate the uncompressed
             * data (size/ CRC) against the gzip footer.
             * Then search for a new stream and, if we
             * find one, re-initialise inflation
             */
            if (z_ret == Z_STREAM_END) {

                zran_log("End of gzip stream [%u]\n", strm->avail_in);

                /*
                 * Make sure we have data in the input buffer
                 * to read and validate the gzip footer and,
                 * in case we are reading concatenated
                 * streams, to search for the next stream and
                 * read its header.
                 *
                 * There is no way of knowing how much data we
                 * need to read in here - there is no upper
                 * bound on the amount of null padding bytes
                 * that may be present in between, or at the
                 * end of, a stream, and there is no upper
                 * bound on the size of a gzip header.

                 * So a critical assumption is made here, that
                 * the size of the read buffer is large enough
                 * to encompass all of:
                 *
                 *   - the footer of a gzip stream,
                 *   - null padding after the end of a gzip
                 *     stream, and
                 *   - the header of the next gzip stream
                 *
                 * This assumption could be removed by changing
                 * the way that data is loaded and buffered
                 * from the file - e.g. we could load bytes in
                 * one-by-one, skipping over null bytes, and
                 * then parse the gzip header ourselves.  But
                 * this is far too much work for what is a very
                 * edge-casey scenario.
                 */
                z_ret = _zran_read_data_from_file(index,
                                                  strm,
                                                  cmp_offset,
                                                  uncmp_offset,
                                                  index->readbuf_size);
                if (!((z_ret == 0) || (z_ret == ZRAN_READ_DATA_EOF))) {
                    goto fail;
                }

                /*
                 * If there is no more data, the footer is
                 * missing, so the data must be corrupt.
                 */
                if (strm->avail_in < 8) {
                    goto fail;
                }

                /*
                 * _validate_stream reads and checks in the
                 * gzip stream footer (the CRC32 and ISIZE
                 * fields at the end of a gzip stream), and
                 * _find_next_stream will skip over any
                 * remaining bytes (e.g. null padding bytes)
                 * until eof, or a new stream is found.
                 *
                 * The number of bytes that were read/
                 * skipped over are accumulated into off.
                 */
                off = 0;

                /*
                 * If we have not yet validated this stream,
                 * check that the CRC and uncompressed size in
                 * the footer match what we have calculated
                 */
                if (uncmp_offset > index->last_stream_ended &&
                    index->validating) {
                    z_ret = _zran_validate_stream(index, strm, &off);

                    if (z_ret == ZRAN_VALIDATE_STREAM_INVALID) {
                        error_return_val = ZRAN_INFLATE_CRC_ERROR;
                        goto fail;
                    }
                    else if (z_ret != 0) {
                        goto fail;
                    }
                    index->last_stream_ended = uncmp_offset;
                    index->validating        = 0;
                }

                /* Otherwise skip over the 8 byte GZIP footer */
                else {
                    strm->avail_in -= 8;
                    strm->next_in  += 8;
                    off            += 8;
                }

                /*
                 * See if we have another concatenated gzip
                 * stream.  If we run out of input data here,
                 * bad things will happen. Refer to the long
                 * comment regarding the input buffer, above.
                 */
                z_ret = _zran_find_next_stream(index, strm, &off);

                cmp_offset      += off;
                _total_consumed += off;

                /*
                 * If _zran_find_next_stream can't find
                 * a new stream, we are either out of
                 * compressed input data, or at eof. In
                 * either case, break and let the outer
                 * loop deal with it.
                 */
                if (z_ret == ZRAN_FIND_STREAM_NOT_FOUND) {
                    break;
                }
                else if (z_ret != 0) {
                    goto fail;
                }

                if (add_stream_points) {
                    if (_zran_add_point(index,
                                        0,
                                        cmp_offset,
                                        uncmp_offset,
                                        0,
                                        0,
                                        NULL) != 0) {
                        goto fail;
                    }
                }
            }

            /*
             * We've run out of space to store decompressed
             * data - this is the responsibility of the caller,
             * so bail out.
             */
            if (strm->avail_out == 0) {

                zran_log("Output buffer full - stopping inflation\n");

                /*
                 * We return OUTPUT_FULL if we haven't
                 * decompressed the requested number of
                 * bytes, or ZRAN_INFLATE_STOP_AT_BLOCK
                 * is active and we haven't yet found a
                 * block.
                 */
                if (inflate_stop_at_block(flags) || _total_output < len) {
                    return_val = ZRAN_INFLATE_OUTPUT_FULL;
                }
                break;
            }
        }

        /*
         * Some of the code above has decided that
         * it wants this _zran_inflate call to return.
         */
        if (return_val != ZRAN_INFLATE_OK) {
            break;
        }
    }

    /*
     * If ZRAN_INFLATE_FREE_READBUF is
     * active, clear input buffer memory
     * and offsets.
     */
    if (inflate_free_readbuf(flags)) {
        free(index->readbuf);
        index->readbuf        = NULL;
        index->readbuf_offset = 0;
        index->readbuf_end    = 0;
    }

    /*
     * Otherwise save the readbuf
     * offset for next time.
     */
    else {
        index->readbuf_offset = index->readbuf_end - strm->avail_in;
    }

    /*
     * If ZRAN_INFLATE_FREE_Z_STREAM
     * is active, do just that.
     */
    if (inflate_free_stream(flags)) {
        if (inflateEnd(strm) != Z_OK)
            goto fail;
    }

    /*
     * Update the total number of
     * bytes that were consumed/read
     */
    *total_consumed = _total_consumed;
    *total_output   = _total_output;

    /*
     * Update the compressed/uncompressed
     * offsets in case we need to use them
     * later.
     */
    index->inflate_cmp_offset   = cmp_offset;
    index->inflate_uncmp_offset = uncmp_offset;

    zran_log("Inflate finished - consumed=%u, output=%u,\n"
             "                   cmp_offset=%llu, uncmp_offset=%llu \n\n",
             *total_consumed, *total_output,
             cmp_offset, uncmp_offset);

    /* Phew. */
    return return_val;

fail:
    if (index->readbuf != NULL) {
        free(index->readbuf);
        index->readbuf        = NULL;
        index->readbuf_offset = 0;
        index->readbuf_end    = 0;
    }

    return error_return_val;
}


/*
 * Expands the index to encompass the
 * compressed offset specified by 'until'.
 */
int _zran_expand_index(zran_index_t *index, uint64_t until) {

    /*
     * Used to store return code when
     * an error occurs.
     */
    int error_return_val = ZRAN_EXPAND_INDEX_FAIL;

    /*
     * Used to store and check return values
     * from zlib and zran functions.
     */
    int z_ret;

    /* Zlib stream struct */
    z_stream strm;

    /*
     * Number of bytes read/decompressed
     * on each call to _zran_inflate.
     */
    uint32_t bytes_consumed;
    uint32_t bytes_output;

    /*
     * Buffer to store uncompressed data,
     * size of said buffer, and current offset
     * into said buffer. We wrap the buffer
     * around to the beginning when it is
     * filled.
     *
     * Ideally, we only want to decompress
     * index->spacing bytes before creating a
     * new index point. But we may have to
     * decompress more than this before a
     * suitable location (a block/stream
     * boundary) is found, so we allocate
     * more space.
     */
    uint8_t *data        = NULL;
    uint32_t data_size   = index->spacing * 4;
    uint32_t data_offset = 0;

    /*
     * _zran_inflate control flags. We need
     * to use different flags on the first
     * call - first_inflate is used to track
     * this.
     */
    uint16_t inflate_flags;
    uint8_t  first_inflate = 1;

    /*
     * Counters to keep track of where we are
     * in both the compressed and uncompressed
     * streams.
     */
    uint64_t cmp_offset;
    uint64_t uncmp_offset;
    uint64_t last_uncmp_offset;

    /*
     * start is a reference to the last
     * point in the index when this function
     * is called. This is where we need
     * to start decompressing data from
     * before we can add more index points.
     *
     * last_created is a reference to the
     * most recent point that was added
     * to the index in this call to
     * _zran_expand_index.
     */
    zran_point_t *start        = NULL;
    zran_point_t *last_created = NULL;

    /*
     * In order to create a new index
     * point, we need to start reading
     * at the last index point, so that
     * we read enough data to initialise
     * the inflation. If we don't have
     * at least two points, we start
     * at the beginning of the file.
     */
    start = NULL;
    if (index->npoints > 1) {

        start = &(index->list[index->npoints - 1]);

        /*
         * The index already covers the requested
         * offset. Nothing needs to be done.
         */
        if (until != 0 && until <= start->cmp_offset)
            return 0;
    }

    /*
     * Allocate memory for the
     * uncompressed data buffer.
     */
    data = calloc(1, data_size);
    if (data == NULL)
        goto fail;

    /* Let's do this. */
    zran_log("_zran_expand_index(%llu)\n", until);

    /*
     * If the caller passed until == 0,
     * we expand to EOF.
     */
    if (until == 0) {
        until = UINT64_MAX;
    }

    /*
     * We start from the last point in
     * the index, or the beginning of
     * the file, if there are not enough
     * points in the index.
     */
    if (start != NULL) {
        cmp_offset        = start->cmp_offset;
        uncmp_offset      = start->uncmp_offset;
        last_uncmp_offset = uncmp_offset;
    }
    else {
        cmp_offset        = 0;
        uncmp_offset      = 0;
        last_uncmp_offset = 0;
    }

    /*
     * Don't finish until we're at the end of the
     * file (break at bottom of loop), or we've
     * expanded the index past the requested
     * offset (and have created at least one new
     * index point - last_created == NULL tells
     * us whether a point has been created).
     */
    while (last_created == NULL || last_created->cmp_offset < until) {

        /*
         * On the first call to _zran_inflate, we
         * tell it to initialise the zlib stream
         * struct, create a read buffer, and start
         * inflation from our starting point.
         */
        if (first_inflate) {
            first_inflate = 0;
            inflate_flags = (ZRAN_INFLATE_INIT_Z_STREAM         |
                             ZRAN_INFLATE_INIT_READBUF          |
                             ZRAN_INFLATE_USE_OFFSET            |
                             ZRAN_INFLATE_CLEAR_READBUF_OFFSETS |
                             ZRAN_INFLATE_STOP_AT_BLOCK);
        }

        /*
         * On subsequent calls, we tell _zran_inflate
         * to just continue where it left off on the
         * previous call.
         */
        else {
            inflate_flags = ZRAN_INFLATE_STOP_AT_BLOCK;
        }

        zran_log("Searching for next block boundary\n"
                 "       c=%llu, u=%llu,\n"
                 "       data_offset=%u, data_space=%u\n",
                 cmp_offset, uncmp_offset, data_offset,
                 data_size - data_offset);

        /*
         * We wrap the data buffer around to its
         * beginning by using some trickery with
         * the data_offset. By doing this, the
         * _zran_add_point function will be able
         * to retrieve the data associated with
         * an index point even if some of it it
         * is contained at the end of the data
         * buffer, and the rest at the beginning.
         */
        z_ret = _zran_inflate(index,
                              &strm,
                              cmp_offset,
                              inflate_flags,
                              &bytes_consumed,
                              &bytes_output,
                              data_size - data_offset,
                              data      + data_offset,
                              1);

        cmp_offset   += bytes_consumed;
        uncmp_offset += bytes_output;
        data_offset   = (data_offset + bytes_output) % data_size;

        /*
         * update the last created offset on every iteration,
         * to catch any index points created by _zran_inflate
         */
        if (index->npoints > 0) {
            last_created      = &index->list[index->npoints - 1];
            last_uncmp_offset = last_created->uncmp_offset;
        }

        /*
         * Has the output buffer been filled?
         * If so, we just continue - the
         * data_offset trickery means that we
         * can ask the _zran_inflate function
         * to just keep filling the buffer
         * until we find a block.
         */
        if (z_ret == ZRAN_INFLATE_OUTPUT_FULL)
            continue;

        /*
         * If z_ret != ZRAN_INFLATE_EOF or
         * ZRAN_INFLATE_BLOCK_BOUNDARY,
         * something has gone wrong.
         */
        else if (z_ret != ZRAN_INFLATE_EOF &&
                 z_ret != ZRAN_INFLATE_BLOCK_BOUNDARY) {
            if (z_ret == ZRAN_INFLATE_CRC_ERROR) {
                error_return_val = ZRAN_EXPAND_INDEX_CRC_ERROR;
            }
            goto fail;
        }

        /*
         * If we're at the end of the file (z_ret
         * == ZRAN_INFLATE_EOF), or at a compress
         * block boundary, and index->spacing bytes
         * have passed since the last index point
         * that was created, we'll create a new
         * index point at this location.
         *
         * Note that the _zran_inflate function
         * also adds index points at the beginning
         * of the file, and at the beginning of all
         * other gzip streams, in the case of
         * concatenated streams (refer to its
         * add_stream_points argument).
         */
        if (z_ret == ZRAN_INFLATE_EOF ||
            uncmp_offset - last_uncmp_offset >= index->spacing) {
            if (_zran_add_point(index,
                                strm.data_type & 7,
                                cmp_offset,
                                uncmp_offset,
                                data_offset,
                                data_size,
                                data) != 0) {
                goto fail;
            }
            last_created      = &index->list[index->npoints - 1];
            last_uncmp_offset = uncmp_offset;
        }

        /* And if at EOF, we are done. */
        if (z_ret == ZRAN_INFLATE_EOF) {
            break;
        }
    }

    /*
     * A final call to _zran_inflate, to clean
     * up read buffer and z_stream memory.
     */
    z_ret = _zran_inflate(index,
                          &strm,
                          0,
                          (ZRAN_INFLATE_CLEAR_READBUF_OFFSETS |
                           ZRAN_INFLATE_FREE_Z_STREAM         |
                           ZRAN_INFLATE_FREE_READBUF),
                          &bytes_consumed,
                          &bytes_output,
                          0,
                          data,
                          0);

    if (z_ret != ZRAN_INFLATE_OK && z_ret != ZRAN_INFLATE_EOF) {
        if (z_ret == ZRAN_INFLATE_CRC_ERROR) {
            error_return_val = ZRAN_EXPAND_INDEX_CRC_ERROR;
        }
        goto fail;
    }

    /*
     * The index may have over-allocated
     * space for storing index points, so
     * here we free the unused memory.
     */
    if (_zran_free_unused(index) != 0) {
        goto fail;
    }

    zran_log("Expansion finished (cmp_offset=%llu, last_created=%llu)\n",
             cmp_offset, last_created->cmp_offset);

    free(data);
    return ZRAN_EXPAND_INDEX_OK;

fail:
    free(data);
    return error_return_val;
}


/*
 * Seek to the approximate location of the specified offset into
 * the uncompressed data stream.
 */
int zran_seek(zran_index_t  *index,
              int64_t        offset,
              uint8_t        whence,
              zran_point_t **point)
{

    int           result;
    zran_point_t *seek_point = NULL;

    zran_log("zran_seek(%lld, %i)\n", offset, whence);

    if (whence == SEEK_END && index->uncompressed_size == 0) {
      goto index_not_built;
    }

    /*
     * The offset passed in is signed, so
     * negative offsets are allowed. But
     * here we transform the offset to
     * positive, as _zran_get_point_with_expand
     * requires an absolute offset from the
     * beginning of the uncompressed stream.
     *
     * I am not currently taking into account
     * the overflow potential when converting
     * from int64 to uint64.
     */

    /*
     * SEEK_END: seek relative to the
     * end of the uncompressed stream
     */
    if (whence == SEEK_END) {
      offset += index->uncompressed_size;
    }

    /*
     * SEEK_CUR: seek relative to
     * the current file position.
     */
    if (whence == SEEK_CUR) {
      offset += index->uncmp_seek_offset;
    }

    /* Bad input */
    if (offset < 0) {
      goto fail;
    }

    /*
     * We implicitly allow seek(0) - if
     * not auto-building the index,
     * seek(0) would otherwwise fail.
     */
    if (offset == 0) {
        index->uncmp_seek_offset = offset;
    }
    else {

        /*
         * Get the index point that
         * corresponds to this offset.
         */
        result = _zran_get_point_with_expand(index, offset, 0, &seek_point);

        if      (result == ZRAN_GET_POINT_EOF)         goto eof;
        else if (result == ZRAN_GET_POINT_NOT_COVERED) goto not_covered;
        else if (result == ZRAN_GET_POINT_CRC_ERROR)   goto crcerror;
        else if (result != ZRAN_GET_POINT_OK)          goto fail;

        /*
         * transform into an offset
         * into the compressed stream
         */
        index->uncmp_seek_offset = offset;
        offset                   = seek_point->cmp_offset;

        /*
         * This index point is not byte-aligned.
         * Adjust the offset accordingly.
         */
        if (seek_point->bits > 0)
            offset -= 1;
    }

    /*
     * The caller wants a ref to the
     * index point corresponding to
     * the seek location.
     */
    if (point != NULL) {
        *point = seek_point;
    }

    if (fseek_(index->fd, index->f, offset, SEEK_SET) != 0)
        goto fail;

    return ZRAN_SEEK_OK;

crcerror:        return ZRAN_SEEK_CRC_ERROR;
fail:            return ZRAN_SEEK_FAIL;
index_not_built: return ZRAN_SEEK_INDEX_NOT_BUILT;
not_covered:     return ZRAN_SEEK_NOT_COVERED;
eof:
    index->uncmp_seek_offset = index->uncompressed_size;
    return ZRAN_SEEK_EOF;
}


/* Return the current seek position in the uncompressed data stream. */
uint64_t zran_tell(zran_index_t *index) {

    return index->uncmp_seek_offset;
}


/* Read len bytes from the uncompressed data stream, storing them in buf. */
int64_t zran_read(zran_index_t *index,
                  void         *buf,
                  uint64_t      len) {

    /* Used to store/check return values. */
    int ret;

    /*
     * Used to store error code for return
     *   if an error occurs
     */
    int error_return_val = ZRAN_READ_FAIL;

    /*
     * Number of bytes we try to read, and
     * number of bytes actually read/output
     * on each call to _zran_inflate.
     */
    uint64_t bytes_to_read;
    uint32_t bytes_consumed;
    uint32_t bytes_output;

    /*
     * _zran_inflate control flags. We need
     * to pass different flags on the first
     * call to _zran_inflate.
     */
    uint16_t inflate_flags;
    uint8_t  first_inflate = 1;

    /*
     * Counters keeping track of the current
     * location in both the compressed and
     * uncompressed streams, and the total
     * number of bytes read.
     */
    uint64_t uncmp_offset;
    uint64_t cmp_offset;
    uint64_t total_read;

    /*
     * Zlib stream struct and starting
     * index point for the read..
     */
    z_stream      strm;
    zran_point_t *start = NULL;

    /*
     * Memory used to store bytes that we skip
     * over before reaching the appropriate
     * point in the uncompressed data stream.
     *
     * to_discard is used to store the number of
     * bytes that we want to discard on a single
     * call to _zran_inflate (which is limited by
     * the discard buffer size).
     *
     * total_discarded keeps track of the total
     * number of bytes discarded so far.
     *
     * discard_size is the size of the discard *
     * buffer. Ideally we will only have to
     * decompress (on average) spacing / 2 bytes
     * before reaching the seek location, but this
     * isn't a guarantee, so we allocate more to
     * reduce the number of reads that are required.
     */
    uint8_t *discard         = NULL;
    uint64_t to_discard      = 0;
    uint64_t total_discarded = 0;
    uint64_t discard_size    = index->spacing * 4;

    if (len == 0)         return 0;
    if (len >  INT64_MAX) goto fail;

    zran_log("zran_read(%llu, %lu)\n", len, index->uncmp_seek_offset);

    /*
     * Search for the index point that corresponds to
     * our current seek location in the uncompressed
     * data stream. Reading from the start of file is
     * always allowed, even if the index does not
     * contain any points.
     */
    if (index->uncmp_seek_offset == 0) {
        cmp_offset   = 0;
        uncmp_offset = 0;
    }
    else {
        ret = _zran_get_point_with_expand(index,
                                          index->uncmp_seek_offset,
                                          0,
                                          &start);

        if      (ret == ZRAN_GET_POINT_EOF)         goto eof;
        if      (ret == ZRAN_GET_POINT_NOT_COVERED) goto not_covered;
        else if (ret != ZRAN_GET_POINT_OK) {
            if (ret == ZRAN_GET_POINT_CRC_ERROR) {
                error_return_val = ZRAN_READ_CRC_ERROR;
            }
            goto fail;
        }

        cmp_offset   = start->cmp_offset;
        uncmp_offset = start->uncmp_offset;
    }

    /*
     * We have to start decompressing from
     * the index point that precedes the seek
     * location, so we need to skip over bytes
     * until we get to that location. We use
     * the discard buffer to store those bytes.
     */
    discard = malloc(discard_size);
    if (discard == NULL) {
        goto fail;
    }

    /*
     * Inflate and discard data until we
     * reach the current seek location
     * into the uncompressed data stream.
     */
    first_inflate   = 1;
    total_discarded = 0;
    while (uncmp_offset < index->uncmp_seek_offset) {

        /*
         * On the first call to _zran_inflate,
         * we tell it to initialise the z_stream,
         * and create a read buffer.
         */
        if (first_inflate) {
            first_inflate = 0;
            inflate_flags = (ZRAN_INFLATE_INIT_Z_STREAM         |
                             ZRAN_INFLATE_INIT_READBUF          |
                             ZRAN_INFLATE_CLEAR_READBUF_OFFSETS |
                             ZRAN_INFLATE_USE_OFFSET);
        }
        /*
         * On subsequent calls, we just tell
         * _zran_inflate to continue where
         * it left off.
         */
        else {
            inflate_flags = 0;
        }

        /*
         * Don't read past the uncompressed seek
         * location - at this point, we will need
         * to stop discarding bytes, and start
         * fulfilling the read request.
         */
        to_discard = index->uncmp_seek_offset - uncmp_offset;
        if (to_discard > discard_size)
            to_discard = discard_size;

        zran_log("Discarding %llu bytes (%llu < %llu)\n",
                 to_discard,
                 uncmp_offset,
                 index->uncmp_seek_offset);

        ret = _zran_inflate(index,
                            &strm,
                            cmp_offset,
                            inflate_flags,
                            &bytes_consumed,
                            &bytes_output,
                            to_discard,
                            discard,
                            0);

        /*
         * _zran_inflate should return 0 if
         * it runs out of output space (which
         * is ok), or it has read enough bytes
         * (which is perfect). Any other
         * return code means that something
         * has gone wrong.
         */
        if (ret != ZRAN_INFLATE_OUTPUT_FULL &&
            ret != ZRAN_INFLATE_EOF         &&
            ret != ZRAN_INFLATE_OK) {
            if (ret == ZRAN_INFLATE_CRC_ERROR) {
                error_return_val = ZRAN_READ_CRC_ERROR;
            }
            goto fail;
        }

        cmp_offset      += bytes_consumed;
        uncmp_offset    += bytes_output;
        total_discarded += bytes_output;
    }

    /*
     * Sanity check - we should be at the
     * correct location in the uncompressed
     * stream.
     *
     * TODO What happens here if we are at EOF?
     */
    if (uncmp_offset != index->uncmp_seek_offset)
        goto fail;

    zran_log("Discarded %llu bytes, ready to "
             "read from %llu (== %llu)\n",
             total_discarded,
             uncmp_offset,
             index->uncmp_seek_offset);

    /*
     * At this point, we are ready to inflate
     * from the uncompressed seek location.
     */

    total_read = 0;
    while (total_read < len) {

        /*
         * If we started at the correct location,
         * the discard loop above will not have
         * executed, and _zran_inflate will not
         * have initialised itself. So we repeat
         * the flag control stuff here.
         */
        if (first_inflate) {
            first_inflate = 0;
            inflate_flags = (ZRAN_INFLATE_INIT_Z_STREAM         |
                             ZRAN_INFLATE_INIT_READBUF          |
                             ZRAN_INFLATE_CLEAR_READBUF_OFFSETS |
                             ZRAN_INFLATE_USE_OFFSET);
        }
        else {
            inflate_flags = 0;
        }

        /*
         * _zran_inflate only allows us to
         * read max(uint32_t) at a time. If
         * len is greater than this, we need
         * to split it into multiple calls.
         */
        bytes_to_read = len - total_read;
        if (bytes_to_read > 4294967295) {
            bytes_to_read = 4294967295;
        }

        ret = _zran_inflate(index,
                            &strm,
                            cmp_offset,
                            inflate_flags,
                            &bytes_consumed,
                            &bytes_output,
                            bytes_to_read,
                            (uint8_t *)(buf) + total_read,
                            0);

        cmp_offset   += bytes_consumed;
        uncmp_offset += bytes_output;
        total_read   += bytes_output;

        if (ret == ZRAN_INFLATE_EOF)
            break;

        else if (ret == ZRAN_INFLATE_OUTPUT_FULL) {

            /*
             * We might be reading 2**32 sized chunks
             * of data on each call to _zran_inflate.
             */
            if (bytes_to_read == len) {
                break;
            }
        }
        else if (ret != ZRAN_INFLATE_OK) {
            if (ret == ZRAN_INFLATE_CRC_ERROR) {
                error_return_val = ZRAN_READ_CRC_ERROR;
            }
            goto fail;
        }

        zran_log("Read %u bytes (%llu / %llu)\n",
                 bytes_output,
                 total_read,
                 len);
    }

    /*
     * A final call to _zran_inflate,
     * to clean up memory
     */
    ret = _zran_inflate(index,
                        &strm,
                        0,
                        (ZRAN_INFLATE_CLEAR_READBUF_OFFSETS |
                         ZRAN_INFLATE_FREE_Z_STREAM         |
                         ZRAN_INFLATE_FREE_READBUF),
                        &bytes_consumed,
                        &bytes_output,
                        0,
                        discard,
                        0);

    if (ret != ZRAN_INFLATE_OK && ret != ZRAN_INFLATE_EOF) {
        if (ret == ZRAN_INFLATE_CRC_ERROR) {
            error_return_val = ZRAN_READ_CRC_ERROR;
        }
        goto fail;
    }

    /*
     * Update the current uncompressed
     * seek position.
     */
    index->uncmp_seek_offset += total_read;

    zran_log("Read succeeded - %llu bytes read [compressed offset: %ld]\n",
             total_read,
             ftell_(index->fd, index->f));

    free(discard);

    return total_read;

not_covered: return ZRAN_READ_NOT_COVERED;
eof:         return ZRAN_READ_EOF;
fail:

    if (discard != NULL)
        free(discard);

    return error_return_val;
}


/*
 * Store checkpoint information from index to file fd. File should be opened
 * in binary write mode.
 */
int zran_export_index(zran_index_t *index,
                      FILE         *fd,
                      PyObject     *f) {

    /*
     * TODO: Endianness check for fwrite calls. Prefer little-endian to be
     * consistent with gzip library.
     */

    /* Used for checking return value of fwrite calls. */
    size_t f_ret;

    /* Used for iterating over elements of zran_index_t.list. */
    zran_point_t *point;
    zran_point_t *list_end;

    /* File flags, currently not used. Also used as a temporary variable. */
    uint8_t flags = 0;

    zran_log("zran_export_index: (%lu, %lu, %u, %u, %u)\n",
             index->compressed_size,
             index->uncompressed_size,
             index->spacing,
             index->window_size,
             index->npoints);

    /* Write ID and version, and check for errors. */
    f_ret = fwrite_(ZRAN_INDEX_FILE_ID, sizeof(ZRAN_INDEX_FILE_ID), 1, fd, f);
    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    f_ret = fwrite_(&ZRAN_INDEX_FILE_VERSION, 1, 1, fd, f);
    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /* Write flags (currently unused) */
    f_ret = fwrite_(&flags, 1, 1, fd, f);
    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /* Write compressed size, and check for errors. */
    f_ret = fwrite_(&index->compressed_size,
                   sizeof(index->compressed_size), 1, fd, f);

    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /* Write uncompressed size, and check for errors. */
    f_ret = fwrite_(&index->uncompressed_size,
                   sizeof(index->uncompressed_size), 1, fd, f);

    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /* Write spacing, and check for errors. */
    f_ret = fwrite_(&index->spacing, sizeof(index->spacing), 1, fd, f);

    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /* Write window size, and check for errors. */
    f_ret = fwrite_(&index->window_size, sizeof(index->window_size), 1, fd, f);

    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /* Write number of points, and check for errors. */
    f_ret = fwrite_(&index->npoints, sizeof(index->npoints), 1, fd, f);

    if (ferror_(fd, f)) goto fail;
    if (f_ret != 1)     goto fail;

    /*
     * We will make two passes over points list now. In the first pass, offset
     * mapping information of each point will be written. In the second pass,
     * checkpoint snapshot data will be written. This will keep offsets bundled
     * together, which enables user to read all offset mappings in one pass.
     */

    /* Write all points iteratively for checkpoint offset mapping. */
    point    = index->list;
    list_end = index->list + index->npoints;
    while (point < list_end) {

        /* Write compressed offset, and check for errors. */
        f_ret = fwrite_(&point->cmp_offset,
                        sizeof(point->cmp_offset), 1, fd, f);
        if (ferror_(fd, f)) goto fail;
        if (f_ret != 1)     goto fail;

        /* Write uncompressed offset, and check for errors. */
        f_ret = fwrite_(&point->uncmp_offset,
                        sizeof(point->uncmp_offset), 1, fd, f);
        if (ferror_(fd, f)) goto fail;
        if (f_ret != 1)     goto fail;

        /* Write bit offset, and check for errors. */
        f_ret = fwrite_(&point->bits, sizeof(point->bits), 1, fd, f);
        if (ferror_(fd, f)) goto fail;
        if (f_ret != 1)     goto fail;

        /* Write data flag, and check for errors. */
        flags = (point->data != NULL) ? 1 : 0;
        f_ret = fwrite_(&flags, 1, 1, fd, f);
        if (ferror_(fd, f)) goto fail;
        if (f_ret != 1)     goto fail;

        zran_log("zran_export_index: (p%lu, %lu, %lu, %u, %u)\n",
                 (index->npoints - (list_end - point)), // point index
                 point->cmp_offset,
                 point->uncmp_offset,
                 point->bits,
                 flags);

        point++;
    }

    /*
     * Now write out the window data for every point. No data is written for
     * points which don't have any data (e.g. at stream boundaries).
     */
    point    = index->list;
    list_end = index->list + index->npoints;
    while (point < list_end) {

        if (point->data == NULL) {
            point++;
            continue;
        }

        /* Write checkpoint data, and check for errors. */
        f_ret = fwrite_(point->data, index->window_size, 1, fd, f);
        if (ferror_(fd, f)) goto fail;
        if (f_ret != 1)     goto fail;

        /* Print first and last three bytes of the checkpoint window. */
        zran_log("zran_export_index: "
                     "(%lu, [%02x %02x %02x...%02x %02x %02x])\n",
                 (index->npoints - (list_end - point)), // point index
                 point->data[0],
                 point->data[1],
                 point->data[2],
                 point->data[index->window_size - 3],
                 point->data[index->window_size - 2],
                 point->data[index->window_size - 1]);

        /* Done with this point. Proceed to next one. */
        point++;
    }

    zran_log("zran_export_index: done\n");

    /*
     * It is important to flush written file when done, since underlying file
     * descriptor can be closed by Python code before having a chance to flush.
     */
    f_ret = fflush_(fd, f);
    if (ferror_(fd, f)) goto fail;
    if (f_ret != 0)     goto fail;

    return ZRAN_EXPORT_OK;

fail:
    return ZRAN_EXPORT_WRITE_ERROR;
}

/*
 * Load checkpoint information from file fd to index. File should be opened in
 * binary read mode.
 */
int zran_import_index(zran_index_t *index,
                      FILE         *fd,
                      PyObject     *f) {

    /* Used for checking return value of fread calls. */
    size_t f_ret;

    /* Return value of function if a failure happens. */
    int fail_ret;

    /* Used for iterating over elements of zran_index_t.list. */
    uint64_t      i;
    zran_point_t *point;
    zran_point_t *list_end;

    /*
     * Used to store flags for each point - allocated once we
     * know how many points there are.
     */
    uint8_t *dataflags = NULL;

    /*
     * Used for checking file ID, version, and
     * flags at the beginning of the file. Flags
     * is also used as a temporary variable.
     */
    char    file_id[sizeof(ZRAN_INDEX_FILE_ID)];
    uint8_t version;
    uint8_t flags;

    /*
     * Data fields that will be read from the file. They aren't stored directly
     * to index struct to keep original index in case of any failures while
     * reading those data.
     */
    uint64_t      compressed_size;
    uint64_t      uncompressed_size;
    uint32_t      spacing;
    uint32_t      window_size;
    uint32_t      npoints;
    zran_point_t *new_list = NULL;

    /* CRC validation is currently not possible on an imported index */
    index->flags |= ZRAN_SKIP_CRC_CHECK;

    /* Check if file is read only. */
    if (!is_readonly(fd, f)) goto fail;

    /* Read ID, and check for file errors and EOF. */
    f_ret = fread_(file_id, sizeof(file_id), 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /* Verify file ID. */
    if (memcmp(file_id, ZRAN_INDEX_FILE_ID, sizeof(file_id)))
        goto unknown_format;

    /* Read file format version */
    f_ret = fread_(&version, 1, 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /* This file is too new for us to cope */
    if (version > ZRAN_INDEX_FILE_VERSION)
        goto unsupported_version;

    /* Read flags (currently unused) */
    f_ret = fread_(&flags, 1, 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /* Read compressed size, and check for file errors and EOF. */
    f_ret = fread_(&compressed_size, sizeof(compressed_size), 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /*
     * Compare compressed_size in the index file to the existing size in
     * the current index (set in zran_init), if they don't match this means
     * this index file is not created for this compressed file.
     */
    if (compressed_size != index->compressed_size)
        goto inconsistent;

    if (feof_(fd, f, f_ret)) goto eof;

    /* Read uncompressed size, and check for file errors and EOF. */
    f_ret = fread_(&uncompressed_size, sizeof(uncompressed_size), 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /*
     * Uncompressed size may not be set in either current index or exported
     * file, or both. Therefore, they are compared only if it's set in both.
     */
    if (uncompressed_size        != 0 &&
        index->uncompressed_size != 0 &&
        index->uncompressed_size != uncompressed_size) goto inconsistent;

    /* Read spacing, and check for file errors and EOF. */
    f_ret = fread_(&spacing, sizeof(spacing), 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /* Read window size, and check for file errors and EOF. */
    f_ret = fread_(&window_size, sizeof(window_size), 1, fd, f);
    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    /*
     * Make sanity checks for window size and spacing. These are similar to
     * sanity checks done in zran_init.
     */
    if (window_size < 32768)       goto fail;
    if (spacing     < window_size) goto fail;

    /* Read number of points, and check for file errors and EOF. */
    f_ret = fread_(&npoints, sizeof(npoints), 1, fd, f);

    if (feof_(fd, f, f_ret)) goto eof;
    if (ferror_(fd, f))      goto read_error;
    if (f_ret != 1)          goto read_error;

    zran_log("zran_import_index: (%u, %lu, %lu, %u, %u, %u)\n",
             version,
             compressed_size,
             uncompressed_size,
             spacing,
             window_size,
             npoints);

    /*
     * At this step, the number of points is known. Allocate space for new list
     * of points. This pointer should be cleaned up before exit in case of
     * failure.
     *
     * The index file is allowed to contain 0 points, in which case we
     * initialise the point list to 8 (same as in zran_init).
     */
    new_list = calloc(1, sizeof(zran_point_t) * max(npoints, 8));
    if (new_list == NULL)
        goto memory_error;

    /*
     * Allocate space for the data flag for each point - whether or not
     * there is data associated with it
     */
    dataflags = calloc(npoints, 1);
    if (dataflags == NULL)
        goto memory_error;

    /* Read new points iteratively for reading offset mapping. */
    for (i = 0, point = new_list; i < npoints; i++, point++) {

        /* Read compressed offset, and check for errors. */
        f_ret = fread_(&point->cmp_offset,
                       sizeof(point->cmp_offset), 1, fd, f);
        if (feof_(fd, f, f_ret)) goto eof;
        if (ferror_(fd, f))      goto read_error;
        if (f_ret != 1)          goto read_error;

        /* Read uncompressed offset, and check for errors. */
        f_ret = fread_(&point->uncmp_offset,
                       sizeof(point->uncmp_offset), 1, fd, f);
        if (feof_(fd, f, f_ret)) goto eof;
        if (ferror_(fd, f))      goto read_error;
        if (f_ret != 1)          goto read_error;

        /* Read bit offset, and check for errors. */
        f_ret = fread_(&point->bits, sizeof(point->bits), 1, fd, f);
        if (feof_(fd, f, f_ret)) goto eof;
        if (ferror_(fd, f))      goto read_error;
        if (f_ret != 1)          goto read_error;

        /* Read data flag (added in version 1), and check for errors. */
        if (version >= 1) {
            f_ret = fread_(&flags, 1, 1, fd, f);
            if (feof_(fd, f, f_ret)) goto eof;
            if (ferror_(fd, f))      goto read_error;
            if (f_ret != 1)          goto read_error;

            /*
             * The data flag determines whether or not any window data
             * is associated with this point. We set point->data to 1
             * to indicate to the loop below that this point has data
             * to be loaded.
             */
        }
        /*
         * In index file version 0, the first point
         * has no data, but all other points do.
         */
        else {
            flags = (point == new_list) ? 0 : 1;
        }

        dataflags[i] = flags;

        zran_log("zran_import_index: (p%lu, %lu, %lu, %u, %u)\n",
                 i,
                 point->cmp_offset,
                 point->uncmp_offset,
                 point->bits,
                 flags);
    }

    /*
     * Now loop through and load the window data for all index points.
     */
    for (i = 0, point = new_list; i < npoints; i++, point++) {

        /*
         * There is no data associated with this point - it is either
         * at the beginning of the file, or on a stream boundary.
         */
        if (dataflags[i] == 0) {
            continue;
        }

        /*
         * Allocate space for checkpoint data. These pointers in each point
         * should be cleaned up in case of any failures.
         */
        point->data = calloc(1, window_size);
        if (point->data == NULL)
            goto memory_error;

        /*
         * Read checkpoint data, and check for errors. End of file can be
         * reached just after the last element, so it's not an error for
         * the last element.
         */
        f_ret = fread_(point->data, window_size, 1, fd, f);
        if (feof_(fd, f, f_ret) && i < npoints - 1) goto eof;
        if (ferror_(fd, f))                         goto read_error;
        if (f_ret != 1)                             goto read_error;

        /*
         * TODO: If there are still more data after importing is done, it
         * is silently ignored. It might be handled by other means.
         */

        /* Print first and last three bytes of the checkpoint window. */
        zran_log("zran_import_index:"
                     "(%lu, [%02x %02x %02x...%02x %02x %02x])\n",
                 i,
                 point->data[0],
                 point->data[1],
                 point->data[2],
                 point->data[window_size - 3],
                 point->data[window_size - 2],
                 point->data[window_size - 1]);
    }

    /* There are no errors, it's safe to overwrite existing index data now. */

    /* If a new uncompressed_size is read, update current index. */
    if (index->uncompressed_size == 0 && uncompressed_size != 0) {
        index->uncompressed_size = uncompressed_size;
    }

    /* Overwrite spacing. */
    if (index->spacing != spacing) {
        index->spacing = spacing;
    }

    /* Overwrite window size. */
    if (index->window_size != window_size) {
        index->window_size = window_size;
    }

    /*
     * Now, we will release current checkpoint list of the index, and then
     * point to the new list.
     */

    /*
     * Initialize point to the second element of the list, and list_end to the
     * end of the list. We initialize to the second element, because first
     * element does not keep checkpoint data.
     */
    point    = index->list + 1;
    list_end = index->list + index->npoints;

    while (point < list_end) {
        free(point->data);
        point++;
    }

    /* Now release the old list. */
    free(index->list);

    /* The old list is dead, long live the new list! */
    index->list    = new_list;
    index->npoints = npoints;

    /*
     * Let's not forget to update the size as well.
     * If npoints is 0, the list will have been
     * initialised to allow space for 8 points.
     */
    index->size    = max(npoints, 8);

    zran_log("zran_import_index: done\n");

    free(dataflags);

    return ZRAN_IMPORT_OK;

    /* For each failure case, we assign return value and then clean up. */
fail:
    fail_ret = ZRAN_IMPORT_FAIL;
    goto cleanup;

eof:
    fail_ret = ZRAN_IMPORT_EOF;
    goto cleanup;

read_error:
    fail_ret = ZRAN_IMPORT_READ_ERROR;
    goto cleanup;

inconsistent:
    fail_ret = ZRAN_IMPORT_INCONSISTENT;
    goto cleanup;

memory_error:
    fail_ret = ZRAN_IMPORT_MEMORY_ERROR;
    goto cleanup;

unknown_format:
    fail_ret = ZRAN_IMPORT_UNKNOWN_FORMAT;
    goto cleanup;

unsupported_version:
    fail_ret = ZRAN_IMPORT_UNSUPPORTED_VERSION;
    goto cleanup;

cleanup:
    if (new_list != NULL) {

        /*
         * Initialize point to the second element of the list, and list_end to
         * the end of the list. We initialize to the second element, because
         * first element does not keep checkpoint data.
         */
        point    = new_list + 1;
        list_end = new_list + npoints;

        /*
         * Release until the end of list or the first NULL data pointer,
         * whichever comes first.
         */
        while (point < list_end && point->data != NULL) {
            free(point->data);
            point++;
        }

        /* Release the list itself. */
        free(new_list);
    }

    if (dataflags != NULL) {
        free(dataflags);
    }

    return fail_ret;
}
