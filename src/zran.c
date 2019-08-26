/*
 * zran.c - indexed access to gzip files.
 *
 * See zran.h for documentation.
 *
 * This module was originally based on the zran example, written by Mark
 * Alder, which ships with the zlib source code.
 *
 * Author: Paul McCarthy <pauldmccarthy@gmail.com>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#ifdef _WIN32
#define FSEEK _fseeki64
#define FTELL _ftelli64
#include "windows.h"
#include "io.h"
static int is_readonly(FILE *fd)
{
    /* Can't find a way to do this correctly under Windows and
       the check is not required anyway since the underlying
       Python module checks it already */
    return 1;
}
#else
#include <fcntl.h>
#define FSEEK fseeko
#define FTELL ftello
/* Check if file is read-only */
static int is_readonly(FILE *fd)
{
    return (fcntl(fileno(fd), F_GETFL) & O_ACCMODE) == O_RDONLY;
}


static uint32_t max(uint32_t a, uint32_t b) {

  if (a > b) return a;
  else       return b;
}
#endif

#include "zran.h"

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

/* Define magic bytes and version for export format. */
const char zran_magic_bytes[] = {'G', 'Z', 'I', 'D', 'X', 0, 0};

/*
 * Discards all points in the index which come after the specfiied
 * compressed offset.
 *
 * Returns 0 on success, non-0 on failure.
 */
static int _zran_invalidate_index(
    zran_index_t *index, /* The index                       */
    uint64_t      from   /* Offset into the compressed data */
);


/*
 * Expands the capacity of the memory used to store the index ilst.
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
  uint8_t       compressed  /* Pass in non-0 to get the compressed stream limit,
                               or 0 for the uncompressed limit. */
);


/* Return codes for _zran_get_point_at */
int ZRAN_GET_POINT_FAIL        =  -1;
int ZRAN_GET_POINT_OK          =  0;
int ZRAN_GET_POINT_NOT_COVERED =  1;
int ZRAN_GET_POINT_EOF         =  2;

/*
 * Searches for the zran_point which preceeds the given offset. The offset
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
 * _zran_get_point_at function, however if the index has been initialised
 * with the ZRAN_AUTO_BUILD flag, this function will never return
 * ZRAN_GET_POINT_NOT_COVERED.
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
 * Used by _zran_expand_index, and _zran_read. Seeks to the specified
 * compressed data offset, and initialises zlib to start
 * decompressing/inflating from said offset.
 *
 * Returns 0 for success, non-0 on failure.
 */
static int _zran_init_zlib_inflate(
    zran_index_t *index,      /* The index */

    z_stream     *stream,     /* Pointer to a z_stream struct */

    zran_point_t *point       /* Pass in NULL to initialise for inflation from
                                 the beginning of the stream. Or pass a
                                 pointer to the index point corresponding to
                                 the location to start from. */
);


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
 * Returns 0 on success, non-0 on failure.
 */
static int _zran_expand_index(
    zran_index_t *index, /* The index                      */
    uint64_t      until  /* Expand the index to this point */
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
                                   data preceeding this index point. */
);


/* _zran_find_next_stream return codes */
int ZRAN_FIND_STREAM_ERROR     = -2;
int ZRAN_FIND_STREAM_NOT_FOUND = -1;

/*
 * This function is used to search for concatenated compressed streams.  It
 * searches through the compressed data (pointed to by stream->next_in) to
 * find the location of the next compressed stream.
 *
 * If a new stream was found, the z_stream struct is re-initialised to
 * decompress data from the new stream. In this case, this function returns
 * the number of bytes in the compressed data that were read before the stream
 * was found.
 *
 * Otherwise (if a compressed stream was not found), this function returns
 * ZRAN_FIND_STREAM_NOT_FOUND.
 *
 * If an error occurs, ZRAN_FIND_STREAM_ERROR is returned.
 */
static int _zran_find_next_stream(
    zran_index_t *index, /* The index           */
    z_stream     *stream /* The z_stream struct */
);


/* _zran_inflate return codes */
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
 * This function may be used in a re-entrant or non-re-entrant manner,
 * depending on the flags which are used. In the latter (more likely) case,
 * various pieces of information representing the current inflation state are
 * stored in fields of the zran_index_t struct.
 *
 * Specifically, this function does the following:
 *
 *   1. Figures out the starting offsets into the compressed/uncompressed
 *      streams. If the ZRAN_INFLATE_USE_OFFSET flag is active, the index
 *      point preceeding the specified offset is used as the starting point.
 *      If there is no such point, ZRAN_INFLATE_NOT_COVERED is returned.  If
 *      ZRAN_INFLATE_USE_OFFSET is not active, index->inflate_cmp_offset and
 *      index->inflate_uncmp_offset are used as the starting point.
 *
 *   2. Initialises the z_stream struct, if
 *      ZRAN_INFLATE_INIT_Z_STREAM is active. Otherwise, the function assumes
 *      that the z_stream struct is already initialised and ready to be used.
 *
 *   3. Create a read buffer, if ZRAN_INFLATE_INIT_READBUF is active. A
 *      reference to the read buffer is stored at index->readbuf.  If
 *      ZRAN_INFLATE_INIT_READBUF is not set, the function assumes that the
 *      read buffer already exists.
 *
 *   4. If the ZRAN_INFLATE_CLEAR_READBUF_OFFSETS flag is active, the read
 *      buffer offset (index->readbuf_offset) and length (index->readbuf_end)
 *      fields are both set to 0. Otherwise, the function assumes that the
 *      current offset/length values are valid.
 *
 *   5. Read some compressed data from the file into the read buffer as needed.
 *
 *   6. Pass that data to the zlib inflate function, and store the resulting
 *      uncompressed data in the provided data buffer.
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
 *   - ZRAN_INFLATE_ERROR:          A critical error has occurred.
 */
static int _zran_inflate(
    zran_index_t *index,          /* Pointer to the index. */

    z_stream     *strm,           /* Pointer to a z_stream struct. */

    uint64_t      offset,         /* Compressed data offset to start inflation
                                     from. */

    uint16_t      flags,          /* Control flags. */

    uint32_t     *total_consumed, /* Pointer which is updated to contain the
                                     total number of bytes that were read
                                     from the input file. */

    uint32_t     *total_output,   /* Pointer which is updated to contain the
                                     total number of bytes that were inflated,
                                     and stored in data. */

    uint32_t      len,            /* Maximum number of bytes to inflate. May
                                     be 0. */

    uint8_t      *data            /* Place to store the inflated bytes. */
);


/* Initialise a zran_index_t struct for use with the given GZIP file. */
int zran_init(zran_index_t *index,
              FILE         *fd,
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
     * The zlib manual specifies that a window size of 32KB is 'always enough'
     * to initialise inflation/deflation with a set dictionary. Less than
     * that is not guaranteed to be enough.
    */
    if (window_size < 32768)
        goto fail;

    /*
     * window_size bytes of uncompressed data are stored with each seek point
     * in the index. So it's a bit silly to have the distance between
     * consecutive points less than the window size.
     */
    if (spacing <= window_size)
      goto fail;

    /* The file must be opened in read-only mode */
    if (!is_readonly(fd))
        goto fail;

    /*
     * Calculate the size of the compressed file
     */
    if (FSEEK(fd, 0, SEEK_END) != 0)
        goto fail;

    compressed_size = FTELL(fd);

    if (compressed_size < 0)
        goto fail;

    if (FSEEK(fd, 0, SEEK_SET) != 0)
        goto fail;

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

    zran_log("_zran_free_unused\n");


    new_list = realloc(index->list, sizeof(zran_point_t) * index->npoints);

    if (new_list == NULL) {
        return -1;
    }

    index->list = new_list;
    index->size = index->npoints;

    return 0;
}


/* Deallocate memory used by a zran_index_t struct. */
void zran_free(zran_index_t *index) {

    uint32_t      i;
    zran_point_t *pt;

    zran_log("zran_free\n");

    for (i = 0; i < index->npoints; i++) {
        pt = &(index->list[i]);

        free(pt->data);
    }

    free(index->list);

    index->fd                = NULL;
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
        return -1;

    if (until == 0)
      until = index->compressed_size;

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
     * uncompressed streams (if the latter is known).
     */
    if (compressed && offset >= index->compressed_size)
        goto eof;

    if (!compressed                  &&
        index->uncompressed_size > 0 &&
        offset >= index->uncompressed_size)
        goto eof;

    zran_log("_zran_get_point_at(%llu, c=%u)\n", offset, compressed);

    /*
     * Figure out how much of the compressed and
     * uncompressed data the index currently covers.
     *
     * No points - no coverage.
     */
    if (index->npoints == 0) {
        cmp_max   = 0;
        uncmp_max = 0;
    }

    /*
     * Otherwise the offsets of the
     * last point in the index.
     */
    else {
        last      = &(index->list[index->npoints - 1]);
        uncmp_max = last->uncmp_offset;
        cmp_max   = last->cmp_offset;
    }

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

    if ((index->flags & ZRAN_AUTO_BUILD) == 0) {
        return _zran_get_point_at(index, offset, compressed, point);
    }

    /*
     * See if there is an index point that
     * covers the specified offset. If there's
     * not, we're going to expand the index
     * until there is.
     */
    result = _zran_get_point_at(index, offset, compressed, point);

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

        /*
         * Expand the index
         */
        if (_zran_expand_index(index, expand) != 0) {
            goto fail;
        }

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
        estimate = round(offset * ((float)last->uncmp_offset / last->cmp_offset));
    }
    else {
        estimate = round(offset * ((float)last->cmp_offset / last->uncmp_offset));
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
     * associated with this point. The
     * first index point (where uncmp_offset == 0)
     * has no data associated with it for
     * obvious reasons.
     */
    if (uncmp_offset == 0) {
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
    if (uncmp_offset > 0) {
        if (data_offset >= index->window_size) {

            memcpy(point_data, data + (data_offset - index->window_size), index->window_size);

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
                            z_stream     *stream,
                            zran_point_t *point)
{

    int     ret;
    int     windowBits;
    int64_t seek_loc;

    windowBits        = index->log_window_size;
    stream->zalloc    = Z_NULL;
    stream->zfree     = Z_NULL;
    stream->opaque    = Z_NULL;
    stream->avail_in  = 0;
    stream->avail_out = 0;
    stream->next_in   = Z_NULL;

    /*
     * Seek to the required location in the compressed
     * data stream. If the provided index point is NULL,
     * we start from the beginning of the file.
     *
     * The compressed offset for index points correspond
     * to the first full byte of compressed data. So if
     * the index point is not byte-aligned (bits > 0), we
     * need to seek to the previous byte, and tell zlib
     * about it (via the inflatePrime call below).
     */
    if (point == NULL) seek_loc = 0;
    else               seek_loc = point->cmp_offset - (point->bits > 0);

    if (FSEEK(index->fd, seek_loc, SEEK_SET) != 0)
        goto fail;

    /*
     * If we're starting from the beginning
     * of the file, we tell inflateInit2 to
     * expect a file header
     */
    if (point == NULL) {

        zran_log("zlib_init_zlib_inflate(0, n/a, n/a, %u + 32)\n", windowBits);
        if (inflateInit2(stream, windowBits + 32) != Z_OK) {
            goto fail;
        }
    }

    /*
     * Otherwise, we configure for raw inflation,
     * and initialise the inflation dictionary
     * from the uncompressed data associated with
     * the index point.
     */
    else {

        zran_log("_zran_init_zlib_inflate(%lld, %llu, %llu, -%u)\n",
                 seek_loc,
                 point->cmp_offset,
                 point->uncmp_offset,
                 windowBits);

        if (inflateInit2(stream, -windowBits) != Z_OK) {
            goto fail;
        }

        /* The starting index point is not
         * byte-aligned, so we'll insert
         * the initial bits into the inflate
         * stream using inflatePrime
         */
        if (point->bits > 0) {

            ret = getc(index->fd);

            if (ret == -1)
                goto fail;

            if (inflatePrime(stream,
                             point->bits, ret >> (8 - point->bits)) != Z_OK)
                goto fail;
        }

        /*
         * Initialise the inflate stream
         * with the index point data.
         */
        if (point->data != NULL) {
            if (inflateSetDictionary(stream,
                                     point->data,
                                     index->window_size) != Z_OK)
                goto fail;
        }
    }

    return 0;

fail:
    return -1;
}


/*
 * Identify the location of the next compressed stream (if the file
 * contains concatenated streams).
 */
int _zran_find_next_stream(zran_index_t *index, z_stream *stream) {

    /*
     * Search for the beginning of
     * the next stream. GZIP files
     * start with 0x1f8b.
     */
    int offset = 0;
    int found  = 0;

    while (stream->avail_in >= 2) {

        if (stream->next_in[0] == 0x1f &&
            stream->next_in[1] == 0x8b) {
            found = 1;
            break;
        }

        offset           += 2;
        stream->next_in  += 2;
        stream->avail_in -= 2;
    }

    /*
     * No header found for
     * the next stream.
     */
    if (found == 0)
        goto not_found;

    zran_log("New stream found, re-initialising inflation\n");

    /*
     * Re-configure for inflation
     * from the new stream.
     */
    if (inflateEnd(stream) != Z_OK)
        goto fail;

    stream->zalloc = Z_NULL;
    stream->zfree  = Z_NULL;
    stream->opaque = Z_NULL;

    if (inflateInit2(stream, index->log_window_size + 32) != Z_OK)
        goto fail;

    return offset;

fail:
    return ZRAN_FIND_STREAM_ERROR;

not_found:
    return ZRAN_FIND_STREAM_NOT_FOUND;
}


/* The workhorse. Inflate/decompress data from the file. */
static int _zran_inflate(zran_index_t *index,
                         z_stream     *strm,
                         uint64_t      offset,
                         uint16_t      flags,
                         uint32_t     *total_consumed,
                         uint32_t     *total_output,
                         uint32_t      len,
                         uint8_t      *data) {

    /*
     * Used to store and check return
     * values. f_ret is for fread,
     * z_ret is for zlib/zran functions.
     * return_val is the return value for
     * this function.
     */
    size_t f_ret;
    int    z_ret;
    int    return_val = ZRAN_INFLATE_OK;

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
     * Index point to start from
     * (if ZRAN_INFLATE_USE_OFFSET
     * is active).
     */
    zran_point_t *start = NULL;

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
     *      point that preceeds the requested offset.
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

    zran_log("initialising to inflate from c=%llu, u=%llu\n",
             cmp_offset,
             uncmp_offset);

    /*
     * If ZRAN_INFLATE_INIT_Z_STREAM is active,
     * initialise the zlib struct for inflation.
     * The _zran_init_zlib_inflate function
     * seeks to the correct location in the file
     * for us.
     *
     * If ZRAN_INFLATE_INIT_Z_STREAM is not
     * active, we assume that the file is
     * already at the correct spot.
     */
    if (inflate_init_stream(flags)) {
        if (_zran_init_zlib_inflate(index, strm, start) != 0) {
            goto fail;
        }
    }

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
     * and index->readbuf_end were sert on a
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
     * Keep going until we run out of space.
     */
    while (strm->avail_out > 0) {

        /*
         * We need to read in more data. We read in more
         * data when strm->avail_in < 2, because a GZIP
         * header is 2 bytes long, and when searching for
         * the next stream in a sequence of concatenated
         * streams, the _zran_find_next_stream function
         * might leave a byte in the input without
         * finding a new stream.
         */
        if (strm->avail_in < 2) {

            /*
             * If there are any unprocessed bytes
             * left over, put them at the beginning
             * of the read buffer
             */
            if (strm->avail_in > 0) {
                memcpy(index->readbuf, strm->next_in, strm->avail_in);
            }

            zran_log("Reading from file %llu [ == %llu?]\n",
                     FTELL(index->fd), cmp_offset);

            /*
             * Read a block of compressed data
             * (offsetting past any left over
             * bytes that we may have copied to
             * the beginning of the read buffer
             * above).
             */
            f_ret = fread(index->readbuf + strm->avail_in,
                          1,
                          index->readbuf_size - strm->avail_in,
                          index->fd);

            if (ferror(index->fd)) {
                goto fail;
            }

            /*
             * No bytes read - we've reached EOF
             */
            if (f_ret == 0) {
                if (feof(index->fd)) {
                    return_val = ZRAN_INFLATE_EOF;
                    break;
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

            zran_log("Read %lu bytes from file [c=%llu, u=%llu]\n",
                     f_ret, cmp_offset, uncmp_offset);

            /*
             * Tell zlib about the block
             * of compressed data that we
             * just read in.
             */
            index->readbuf_end = f_ret + strm->avail_in;
            strm->avail_in    += f_ret;
            strm->next_in      = index->readbuf;
        }

        /*
         * Decompress the block until it is
         * gone (or we've read enough bytes)
         */
        z_ret = Z_OK;
        while (strm->avail_in > 0) {

            /*
             * Re-initialise inflation if we have
             * hit a new compressed stream.
             */
            if (z_ret == Z_STREAM_END) {

                zran_log("End of stream - searching for another stream\n");

                z_ret = _zran_find_next_stream(index, strm);

                /*
                 * If _zran_find_next_stream can't find
                 * a new stream, we are either out of
                 * compressed input data, or at eof. In
                 * either case, break and let the outer
                 * loop deal with it.
                 */
                if      (z_ret == ZRAN_FIND_STREAM_NOT_FOUND) break;
                else if (z_ret == ZRAN_FIND_STREAM_ERROR)     goto fail;

                /*
                 * _zran_find_next_stream has found a
                 * new stream, and has told us how many
                 * bytes it skipped over before finding
                 * it.
                 */
                cmp_offset      += z_ret;
                _total_consumed += z_ret;
            }

            /*
             * Optimistically update offsets -
             * we will adjust them after the
             * inflate call.
             */
            cmp_offset      += strm->avail_in;
            uncmp_offset    += strm->avail_out;
            _total_consumed += strm->avail_in;
            _total_output   += strm->avail_out;

            zran_log("Before inflate - avail_in=%u, avail_out=%u\n",
                     strm->avail_in, strm->avail_out);

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
            if (inflate_stop_at_block(flags)) z_ret = inflate(strm, Z_BLOCK);
            else                              z_ret = inflate(strm, Z_NO_FLUSH);

            zran_log("After inflate - avail_in=%u, avail_out=%u\n",
                     strm->avail_in, strm->avail_out);

            /*
             * Adjust our offsets according to what
             * was actually consumed/decompressed.
             */
            cmp_offset      -= strm->avail_in;
            uncmp_offset    -= strm->avail_out;
            _total_consumed -= strm->avail_in;
            _total_output   -= strm->avail_out;

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
             * We've run out of space to
             * store decompressed data
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

            /*
             * End of file. The GZIP file
             * footer takes up 8 bytes, which
             * do not get processed by the
             * inflate function.
             *
             * We use ftell rather than feof,
             * as the EOF indicator only gets
             * set on attempts to read past
             * the end of a file, and this
             * won't happen when the file
             * size is an exact multiple of
             # the read buffer size.
             */
            if ((FTELL(index->fd) >= index->compressed_size) &&
                strm->avail_in <= 8) {

                zran_log("End of file, stopping inflation\n");

                return_val = ZRAN_INFLATE_EOF;

                /*
                 * We now know how big the
                 * uncompressed data is.
                 */
                if (index->uncompressed_size == 0) {

                    zran_log("Updating uncompressed data "
                             "size: %llu\n", uncmp_offset);
                    index->uncompressed_size = uncmp_offset;
                }
                break;
            }

            /*
             * Some of the code above has decided that
             * it wants this _zran_inflate call to return.
             */
            if (return_val != ZRAN_INFLATE_OK) {
                break;
            }
        }

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

    return ZRAN_INFLATE_ERROR;
}


/*
 * Expands the index to encompass the
 * compressed offset specified by 'until'.
 */
int _zran_expand_index(zran_index_t *index, uint64_t until)
{

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
     * into said buffef. We wrap the buffer
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
        if (until <= start->cmp_offset)
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
     * we force some data to be read.
     */
    if (until == 0) {
      until = index->spacing;
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
     * file, or we've expanded the index past
     * the requested offset (and have created at
     * least one new index point -
     * last_created == NULL tells us whether a
     * point has been created).
     */
    while ((cmp_offset < index->compressed_size) &&
           (last_created == NULL || last_created->cmp_offset < until)) {

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
                              data      + data_offset);

        cmp_offset   += bytes_consumed;
        uncmp_offset += bytes_output;
        data_offset   = (data_offset + bytes_output) % data_size;

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
            goto fail;
        }

        /*
         * If we're at the beginning of the file
         * (uncmp_offset == 0), or at the end of
         * the file (z_ret == ZRAN_INFLATE_EOF),
         * or at a compress block boundary,
         * and index->spacing bytes have passed
         * since the last index point that was
         * created, we'll create a  new index
         * point at this location.
         */
        if (z_ret == ZRAN_INFLATE_EOF ||
            uncmp_offset == 0         ||
            uncmp_offset - last_uncmp_offset >= index->spacing) {

            // TODO If at start or EOF, you should
            //      pass in  NULL for the window
            //      data. You can then clean up
            //      _zran_add_point a little bit.

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
        if (z_ret == ZRAN_INFLATE_EOF)
            break;
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
                          data);

    if (z_ret != ZRAN_INFLATE_OK && z_ret != ZRAN_INFLATE_EOF) {
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
    return 0;

fail:
    free(data);

    return -1;
}


/*
 * Seek to the approximate location of the specified offset into
 * the uncompressed data stream. The whence argument must be
 * SEEK_SET or SEEK_CUR.
 */
int zran_seek(zran_index_t  *index,
              int64_t        offset,
              uint8_t        whence,
              zran_point_t **point)
{

    int           result;
    zran_point_t *seek_point;

    zran_log("zran_seek(%lld, %i)\n", offset, whence);

    if (whence != SEEK_SET && whence != SEEK_CUR) {
        goto fail;
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
     * Get the index point that
     * corresponds to this offset.
     */
    result = _zran_get_point_with_expand(index, offset, 0, &seek_point);

    if (result == ZRAN_GET_POINT_FAIL)        goto fail;
    if (result == ZRAN_GET_POINT_NOT_COVERED) goto not_covered;
    if (result == ZRAN_GET_POINT_EOF)         goto eof;

    index->uncmp_seek_offset = offset;
    offset                   = seek_point->cmp_offset;

    /*
     * This index point is not byte-aligned.
     * Adjust the offset accordingly.
     */
    if (seek_point->bits > 0)
        offset -= 1;

    /*
     * The caller wants a ref to the
     * index point corresponding to
     * the seek location.
     */
    if (point != NULL) {
        *point = seek_point;
    }

    if (FSEEK(index->fd, offset, SEEK_SET) != 0)
        goto fail;

    return ZRAN_SEEK_OK;

fail:        return ZRAN_SEEK_FAIL;
not_covered: return ZRAN_SEEK_NOT_COVERED;
eof:
    index->uncmp_seek_offset = index->uncompressed_size;
    return ZRAN_SEEK_EOF;
}

/* Return the current seek position in the uncompressed data stream. */
long zran_tell(zran_index_t *index) {

    return index->uncmp_seek_offset;
}


/* Read len bytes from the uncompressed data stream, storing them in buf. */
int64_t zran_read(zran_index_t *index,
                  void         *buf,
                  uint64_t      len) {

    /* Used to store/check return values. */
    int ret;

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
     * to pass different flags on thefirst
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

    zran_log("zran_read(%llu)\n", len);

    /*
     * Search for the index point that
     * corresponds to our current seek
     * location in the uncompressed
     * data stream.
     */
    ret = _zran_get_point_with_expand(index,
                                      index->uncmp_seek_offset,
                                      0,
                                      &start);

    if (ret == ZRAN_GET_POINT_EOF)         goto eof;
    if (ret == ZRAN_GET_POINT_NOT_COVERED) goto not_covered;

    /*
     * We have to start decompressing from
     * the index point that preceeds the seek
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
     * into the uncompresesd data stream.
     */
    cmp_offset      = start->cmp_offset;
    uncmp_offset    = start->uncmp_offset;
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
                            discard);

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
            ret != ZRAN_INFLATE_OK)
            goto fail;

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
                            (uint8_t *)(buf) + total_read);

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
        else if (ret != ZRAN_INFLATE_OK)
            goto fail;

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
                        discard);

    if (ret != ZRAN_INFLATE_OK && ret != ZRAN_INFLATE_EOF) {
        goto fail;
    }

    /*
     * Update the current uncompressed
     * seek position.
     */
    index->uncmp_seek_offset += total_read;

    zran_log("Read succeeded - %llu bytes read [compressed offset: %ld]\n",
             total_read,
             FTELL(index->fd));

    free(discard);

    return total_read;

not_covered: return ZRAN_READ_NOT_COVERED;
eof:         return ZRAN_READ_EOF;

fail:

    if (discard != NULL)
        free(discard);

    return ZRAN_READ_FAIL;
}

/*
 * Store checkpoint information from index to file fd. File should be opened in
 * binary write mode.
 */
int zran_export_index(zran_index_t *index,
                      FILE *fd) {

    /*
     * TODO: Endianness check for fwrite calls. Prefer little-endian to be
     * consistent with gzip library.
     */

    /* Used for checking return value of fwrite calls. */
    size_t f_ret;

    /* Used for iterating over elements of zran_index_t.list. */
    zran_point_t *point;
    zran_point_t *list_end;

    /*
     * compressed_size and uncompressed_size are defined as size_t in
     * zran_index_t. They are stored in fixed-length 64-bit variables to be
     * exported portably. Other fields in zran_index_t are defined with
     * fixed-length types, so they are exported directly from index.
     */
    uint64_t compressed_size   = index->compressed_size;
    uint64_t uncompressed_size = index->uncompressed_size;

    zran_log("zran_export_index: (%lu, %lu, %u, %u, %u)\n",
             index->compressed_size,
             index->uncompressed_size,
             index->spacing,
             index->window_size,
             index->npoints);

    /* Write magic bytes, and check for errors. */
    f_ret = fwrite(zran_magic_bytes, sizeof(zran_magic_bytes), 1, fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 1) goto fail;

    /* Write compressed size, and check for errors. */
    f_ret = fwrite(&compressed_size, sizeof(compressed_size), 1, fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 1) goto fail;

    /* Write uncompressed size, and check for errors. */
    f_ret = fwrite(&uncompressed_size, sizeof(uncompressed_size), 1, fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 1) goto fail;

    /* Write spacing, and check for errors. */
    f_ret = fwrite(&index->spacing, sizeof(index->spacing), 1, fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 1) goto fail;

    /* Write window size, and check for errors. */
    f_ret = fwrite(&index->window_size, sizeof(index->window_size), 1, fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 1) goto fail;

    /* Write number of points, and check for errors. */
    f_ret = fwrite(&index->npoints, sizeof(index->npoints), 1, fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 1) goto fail;

    /*
     * We will make two passes over points list now. In the first pass, offset
     * mapping information of each point will be written. In the second pass,
     * checkpoint snapshot data will be written. This will keep offsets bundled
     * together, which enables user to read all offset mappings in one pass.
     */

    /*
     * Initialize point to the first element of the list, and list_end to the
     * end of the list.
     */
    point    = index->list;
    list_end = index->list + index->npoints;

    /* Write all points iteratively for checkpoint offset mapping. */
    while (point < list_end) {

        /*
         * Below, it's preferred to write members of points elementwise to
         * avoid handling struct alignment and padding issues.
         */

        /* Write compressed offset, and check for errors. */
        f_ret = fwrite(&point->cmp_offset, sizeof(point->cmp_offset), 1, fd);

        if (ferror(fd)) goto fail;
        if (f_ret != 1) goto fail;

        /* Write uncompressed offset, and check for errors. */
        f_ret = fwrite(&point->uncmp_offset, sizeof(point->uncmp_offset), 1, fd);

        if (ferror(fd)) goto fail;
        if (f_ret != 1) goto fail;

        /* Write bit offset, and check for errors. */
        f_ret = fwrite(&point->bits, sizeof(point->bits), 1, fd);

        if (ferror(fd)) goto fail;
        if (f_ret != 1) goto fail;

        zran_log("zran_export_index: (%lu, %lu, %lu, %u)\n",
                 (index->npoints - (list_end - point)), // point index
                 point->cmp_offset,
                 point->uncmp_offset,
                 point->bits);

        /* Done with this point. Proceed to next one. */
        point++;
    }

    /*
     * Initialize point to the second element of the list (as first element has
     * no data), and list_end to the end of the list.
     */
    point    = index->list + 1;
    list_end = index->list + index->npoints;

    /* Write all points iteratively for writing checkpoint data. */
    while (point < list_end) {

        /* Write checkpoint data, and check for errors. */
        f_ret = fwrite(point->data, index->window_size, 1, fd);

        if (ferror(fd)) goto fail;
        if (f_ret != 1) goto fail;

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
    f_ret = fflush(fd);

    if (ferror(fd)) goto fail;
    if (f_ret != 0) goto fail;

    return ZRAN_EXPORT_OK;

fail:
    return ZRAN_EXPORT_WRITE_ERROR;

}

/*
 * Load checkpoint information from file fd to index. File should be opened in
 * binary read mode.
 */
int zran_import_index(zran_index_t *index,
                      FILE *fd) {

    /* Used for checking return value of fread calls. */
    size_t f_ret;

    /* Return value of function if a failure happens. */
    int fail_ret;

    /* Used for iterating over elements of zran_index_t.list. */
    zran_point_t *point;
    zran_point_t *list_end;

    /* Used for checking magic bytes in the beginning of the file. */
    char magic_bytes[sizeof(zran_magic_bytes)];

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

    /* Check if file is read only. */
    if (!is_readonly(fd)) goto fail;

    /* Read magic bytes, and check for file errors and EOF. */
    f_ret = fread(magic_bytes, sizeof(magic_bytes), 1, fd);

    if (feof(fd))   goto eof;
    if (ferror(fd)) goto read_error;
    if (f_ret != 1) goto read_error;

    /* Verify magic bytes. */
    if (memcmp(magic_bytes, zran_magic_bytes, sizeof(magic_bytes))) {
        goto unknown_format;
    }

    /* Read compressed size, and check for file errors and EOF. */
    f_ret = fread(&compressed_size, sizeof(compressed_size), 1, fd);

    if (feof(fd))   goto eof;
    if (ferror(fd)) goto read_error;
    if (f_ret != 1) goto read_error;

    /*
     * Make sanity checks for compressed size. Since compressed_size is of type
     * size_t, and the value encoded in the export file is 64-bit, it can
     * overflow. It is also compared to the existing size in the current index
     * (set in zran_init), if they don't match this means this index file is
     * not created for this compressed file.
     */
    if (compressed_size >= SIZE_MAX)               goto overflow;
    if (compressed_size != index->compressed_size) goto inconsistent;

    /* Read uncompressed size, and check for file errors and EOF. */
    f_ret = fread(&uncompressed_size, sizeof(uncompressed_size), 1, fd);

    if (feof(fd))   goto eof;
    if (ferror(fd)) goto read_error;
    if (f_ret != 1) goto read_error;

    /*
     * Make sanity checks for uncompressed size. Similar to compressed_size,
     * the value is checked against overflow. Uncompressed size may not be set
     * in either current index or exported file, or both. Therefore, they are
     * compared only if it's set in both.
     */
    if (uncompressed_size >= SIZE_MAX) goto overflow;
    if (uncompressed_size        != 0 &&
        index->uncompressed_size != 0 &&
        index->uncompressed_size != uncompressed_size) goto inconsistent;

    /* Read spacing, and check for file errors and EOF. */
    f_ret = fread(&spacing, sizeof(spacing), 1, fd);

    if (feof(fd))   goto eof;
    if (ferror(fd)) goto read_error;
    if (f_ret != 1) goto read_error;

    /* Read window size, and check for file errors and EOF. */
    f_ret = fread(&window_size, sizeof(window_size), 1, fd);

    if (feof(fd))   goto eof;
    if (ferror(fd)) goto read_error;
    if (f_ret != 1) goto read_error;

    /*
     * Make sanity checks for window size and spacing. These are similar to
     * sanity checks done in zran_init.
     */
    if (window_size < 32768)       goto fail;
    if (spacing     < window_size) goto fail;

    /* Read number of points, and check for file errors and EOF. */
    f_ret = fread(&npoints, sizeof(npoints), 1, fd);

    if (feof(fd))   goto eof;
    if (ferror(fd)) goto read_error;
    if (f_ret != 1) goto read_error;

    zran_log("zran_import_index: (%lu, %lu, %u, %u, %u)\n",
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

    if (new_list == NULL) goto memory_error;

    /*
     * Initialize point to the first element of the list, and list_end to the
     * end of the list.
     */
    point    = new_list;
    list_end = new_list + npoints;

    /* Read new points iteratively for reading offset mapping. */
    while (point < list_end)
    {

        /* Read compressed offset, and check for errors. */
        f_ret = fread(&point->cmp_offset, sizeof(point->cmp_offset), 1, fd);

        if (feof(fd))   goto eof;
        if (ferror(fd)) goto read_error;
        if (f_ret != 1) goto read_error;

        /* Read uncompressed offset, and check for errors. */
        f_ret = fread(&point->uncmp_offset, sizeof(point->uncmp_offset), 1, fd);

        if (feof(fd))   goto eof;
        if (ferror(fd)) goto read_error;
        if (f_ret != 1) goto read_error;

        /* Read bit offset, and check for errors. */
        f_ret = fread(&point->bits, sizeof(point->bits), 1, fd);

        if (feof(fd))   goto eof;
        if (ferror(fd)) goto read_error;
        if (f_ret != 1) goto read_error;

        zran_log("zran_import_index: (%lu, %lu, %lu, %u)\n",
                 (npoints - (list_end - point)), // point index
                 point->cmp_offset,
                 point->uncmp_offset,
                 point->bits);

        /* Done with this point. Proceed to the next one. */
        point++;
    }

    /*
     * Initialize point to the second element of the list (as the first element
     * has no data), and list_end to the end of the list.
     */
    point    = new_list + 1;
    list_end = new_list + npoints;

    /* Read new points iteratively for reading checkpoint data. */
    while (point < list_end) {

        /*
         * Allocate space for checkpoint data. These pointers in each point
         * should be cleaned up in case of any failures.
         */
        point->data = calloc(1, window_size);

        if (point->data == NULL) goto memory_error;

        /*
         * Read checkpoint data, and check for errors. End of file can be
         * reached just after the last element, so it's not an error for
         * the last element.
         */
        f_ret = fread(point->data, window_size, 1, fd);

        if (feof(fd) && point < list_end - 1) goto eof;
        if (ferror(fd))                       goto read_error;
        if (f_ret != 1)                       goto read_error;

        /*
         * TODO: If there are still more data after importing is done, it
         * is silently ignored. It might be handled by other means.
         */

        /* Print first and last three bytes of the checkpoint window. */
        zran_log("zran_import_index:"
                     "(%lu, [%02x %02x %02x...%02x %02x %02x])\n",
                 (npoints - (list_end - point)), // point index
                 point->data[0],
                 point->data[1],
                 point->data[2],
                 point->data[window_size - 3],
                 point->data[window_size - 2],
                 point->data[window_size - 1]);

        /* Done with this point. Proceed to the next one. */
        point++;
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

overflow:
    fail_ret = ZRAN_IMPORT_OVERFLOW;
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

    return fail_ret;
}
