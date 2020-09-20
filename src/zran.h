#ifndef __ZRAN_H__
#define __ZRAN_H__

/*
 * The zran module is an adaptation of the zran example, written by Mark
 * Alder, which ships with the zlib source code. It allows the creation
 * of an index into a compressed file, which is used to improve the speed
 * of random seek/read access to the uncompressed data.
 *
 * Author: Paul McCarthy <pauldmccarthy@gmail.com>
 */

#include <stdlib.h>
#include <stdint.h>


struct _zran_index;
struct _zran_point;


typedef struct _zran_index zran_index_t;
typedef struct _zran_point zran_point_t;

/*
 * These values may be passed in as flags to the zran_init function.
 * They are specified as bit-masks, rather than bit locations.
 */
enum {
  ZRAN_AUTO_BUILD = 1,

};


/*
 * Struct representing the index. None of the fields in this struct
 * should ever need to be accessed or modified directly.
 */
struct _zran_index {

    /*
     * Handle to the compressed file.
     */
    FILE         *fd;

    /*
     * Size of the compressed file. This
     * is calculated in zran_init.
     */
    uint64_t      compressed_size;

    /*
     * Size of the uncompressed data. This is
     * only updated when it becomes known.
     */
    uint64_t      uncompressed_size;

    /*
     * Spacing size in bytes, relative to the
     * uncompressed data stream, between adjacent
     * index points.
     */
    uint32_t      spacing;

    /*
     * Number of bytes of uncompressed data to store
     * for each index point. This must be a minimum
     * of 32768 bytes.
     */
    uint32_t      window_size;

    /*
     * Base2 logarithm of the window size - it
     * is needed to initialise zlib inflation.
     */
    uint32_t      log_window_size;

    /*
     * Size, in bytes, of buffer used to store
     * compressed data read from disk.
     */
    uint32_t      readbuf_size;

    /*
     * Number of index points that have been created.
     */
    uint32_t      npoints;

    /*
     * Number of index points that can be stored.
     */
    uint32_t      size;

    /*
     * List of index points.
     */
    zran_point_t *list;

    /*
     * Most recently requested seek/read
     * location into the uncompressed data
     * stream - this is used to keep track
     * of where the calling code thinks it
     * is in the (uncompressed) file.
     */
    uint64_t      uncmp_seek_offset;

    /*
     * Flags passed to zran_init
     */
    uint16_t      flags;

    /*
     * All of the fields after this point are used
     * by the internal _zran_inflate function.
     */

    /*
     * Reference to a file input
     * buffer of size readbuf_size.
     */
    uint8_t      *readbuf;

    /*
     * An offset into readbuf.
     */
    uint32_t      readbuf_offset;

    /*
     * The current end of the readbuf contents.
     */
    uint32_t      readbuf_end;

    /*
     * Current offsets into the uncompressed and
     * compressed data streams.
     */
    uint64_t      inflate_cmp_offset;
    uint64_t      inflate_uncmp_offset;
};


/*
 * Struct representing a single seek point in the index.
 */
struct _zran_point {


    /*
     * Location of this point in the compressed data
     * stream. This is the location of the first full
     * byte of compressed data - if  the compressed
     * and uncompressed locations are not byte-aligned,
     * the bits field below specifies the bit offset.
     */
    uint64_t  cmp_offset;

    /*
     * Corresponding location of this point
     * in the uncompressed data stream.
     */
    uint64_t  uncmp_offset;

    /*
     * If this point is not byte-aligned, this specifies
     * the number of bits, in the compressed stream,
     * back from cmp_offset, that the uncompressed data
     * starts.
     */
    uint8_t   bits;

    /*
     * Chunk of uncompressed data preceeding this point.
     * This is required to initialise decompression from
     * this point onwards.
     */
    uint8_t  *data;
};


/*
 * Initialise a zran_index_t struct for use with the given file.
 *
 * Passing in 0 for the spacing, window_size and readbuf_size arguments
 * will result in the follwoing values being used:
 *
 *    spacing:      1048576
 *    window_size:  32768
 *    readbuf_size: 16384
 *
 * The flags argument is a bit mask used to control the following options:
 *
 *     ZRAN_AUTO_BUILD: Build the index automatically on demand.
 */
int  zran_init(
  zran_index_t *index,        /* The index                          */
  FILE         *fd,           /* Open handle to the compressed file */
  uint32_t      spacing,      /* Distance in bytes between
                                 index seek points                  */
  uint32_t      window_size,  /* Number of uncompressed bytes
                                 to store with each point           */
  uint32_t      readbuf_size, /* Number of bytes to read at a time  */
  uint16_t      flags         /* Flags controlling index behaviour  */
);


/*
 * Frees the memory use by the given index. The zran_index_t struct
 * itself is not freed.
 */
void zran_free(
  zran_index_t *index /* The index */
);


/*
 * (Re-)Builds the index to cover the given range, which must be
 * specified relative to the compressed data stream. Pass in 0
 * for both offsets to re-build the full index.
 *
 * Returns 0 on success, non-0 on failure.
 */
int zran_build_index(
  zran_index_t *index, /* The index */
  uint64_t      from,  /* Build the index from this point */
  uint64_t      until  /* Build the index to this point   */
);


/* Return codes for zran_seek. */
enum {
    ZRAN_SEEK_FAIL            = -1,
    ZRAN_SEEK_OK              =  0,
    ZRAN_SEEK_NOT_COVERED     =  1,
    ZRAN_SEEK_EOF             =  2,
    ZRAN_SEEK_INDEX_NOT_BUILT =  3
};


/*
 * Seek to the specified offset in the uncompressed data stream.
 * If the index does not currently cover the offset, and it was
 * created with the ZRAN_AUTO_BUILD flag, the index is expanded
 * to cover the offset.
 *
 * Seeking from the end of the uncompressed stream (using SEEK_END)
 * is only possible if the index fully covers the file.
 *
 * Returns:
 *    - ZRAN_SEEK_OK for success.
 *
 *    - ZRAN_SEEK_INDEX_NOT_BUILT if SEEK_END is used, and the index
 *      does not fully cover the file.
 *
 *    - ZRAN_SEEK_NOT_COVERED to indicate that the index does not
 *      cover the requested offset (will never happen if
 *      ZRAN_AUTO_BUILD is active).
 *
 *    - ZRAN_SEEK_EOF to indicate that the requested offset
 *      is past the end of the uncompressed stream.
 *
 *    - ZRAN_SEEK_FAIL to indicate failure of some sort.
 */
int zran_seek(
  zran_index_t  *index,   /* The index                       */
  int64_t        offset,  /* Uncompressed offset to seek to  */
  uint8_t        whence,  /* SEEK_SET, SEEK_CUR, or SEEK_END */
  zran_point_t **point    /* Optional place to store
                             corresponding zran_point_t      */
);

/*
 * Returns the current seek location in the uncompressed data stream
 * (just returns zran_index_t.uncmp_seek_offset).
 */
uint64_t zran_tell(
  zran_index_t *index /* The index */
);


/* Return codes for zran_read. */
enum {
    ZRAN_READ_NOT_COVERED = -1,
    ZRAN_READ_EOF         = -2,
    ZRAN_READ_FAIL        = -3
};

/*
 * Read len bytes from the current location in the uncompressed
 * data stream, storing them in buf. If the index was created with
 * the ZRAN_AUTO_BUILD flag, it is expanded as needed.
 *
 * Returns:
 *   - Number of bytes read for success.
 *
 *   - ZRAN_READ_NOT_COVERED to indicate that the index does not
 *     cover the requested region (will never happen if
 *     ZRAN_AUTO_BUILD is active).
 *
 *   - ZRAN_READ_EOF to indicate that the read could not be completed
 *     because the current uncompressed seek point is at EOF.
 *
 *   - ZRAN_READ_FAIL to indicate that the read failed for some reason.
 */
int64_t zran_read(
  zran_index_t  *index, /* The index                 */
  void          *buf,   /* Buffer to store len bytes */
  uint64_t       len    /* Number of bytes to read   */
);

/* Magic bytes for exported index type. Value is defined in zran.c. */
extern const char zran_magic_bytes[];

/* Return codes for zran_export_index. */
enum {
    ZRAN_EXPORT_OK          =  0,
    ZRAN_EXPORT_WRITE_ERROR = -1
};

/*
 * Export current index data to given file. This exported file later can be
 * used to rebuild index without needing to going through the file again.
 *
 * See zran_import_index for importing.
 *
 * A zran index file is a binary file which has the following header
 * structure. All fields are assumed to be stored with little-endian
 * ordering:
 *
 * | Offset | Length | Description                     |
 * | 0      | 7      | File header (GZIDX\00\00)       |
 * | 7      | 8      | Compressed file size  (uint64)  |
 * | 15     | 8      | Uncompressed file size (uint64) |
 * | 23     | 4      | Index point spacing (uint32)    |
 * | 27     | 4      | Index window size W (uint32)    |
 * | 31     | 4      | Number of index points (uint32) |
 *
 * The header is followed by the offsets for each index point:
 *
 * | Offset | Length | Description                              |
 * | 0      | 8      | Compressed offset for point 0 (uint64)   |
 * | 8      | 8      | Uncompressed offset for point 0 (uint64) |
 * | 16     | 1      | Bit offset for point 0 (uint8)           |
 * | ...    | ...    | ...                                      |
 * | N*17   | 8      | Compressed offset for point N (uint64)   |
 * | ...    | ...    | ...                                      |
 *
 * Finally the window data for every index point is concatenated
 * (W represents the index window size):
 *
 * | Offset | Length | Description                   |
 * | 0      | W      | Window data for index point N |
 * | ...    | ...    | ...                           |
 * | N*W    | W      | Window data for index point N |
 *
 * Returns:
 *   - ZRAN_EXPORT_OK for success.
 *
 *   - ZRAN_EXPORT_WRITE_ERROR to indicate an error from writing to underlying
 *     file.
 */
int zran_export_index(
  zran_index_t  *index, /* The index                  */
  FILE          *fd     /* Open handle to export file */
);

/* Return codes for zran_import_index. */
enum {
    ZRAN_IMPORT_OK             =  0,
    ZRAN_IMPORT_FAIL           = -1,
    ZRAN_IMPORT_EOF            = -2,
    ZRAN_IMPORT_READ_ERROR     = -3,
    ZRAN_IMPORT_INCONSISTENT   = -4,
    ZRAN_IMPORT_MEMORY_ERROR   = -5,
    ZRAN_IMPORT_UNKNOWN_FORMAT = -6
};

/*
 * Import current index from the given file.  index must have been initialized
 * by zran_init function before calling this function, as it is not supported
 * importing into an uninitialised zran_index_t struct. Existing index will be
 * overwritten including spacing and window_size values, whereas values of
 * readbuf_size and flags will be kept.
 *
 * Updating an index file is not supported currently. To update an index file,
 * first import it, create new checkpoints, and then export it again.
 *
 * See zran_export_index for exporting.
 *
 * Returns:
 *   - ZRAN_IMPORT_OK for success.
 *
 *   - ZRAN_IMPORT_FAIL general errors.
 *
 *   - ZRAN_IMPORT_EOF to indicate unexpected end-of-file.
 *
 *   - ZRAN_IMPORT_READ_ERROR to indicate error while reading file.
 *
 *   - ZRAN_IMPORT_OVERFLOW to indicate overflow while reading compressed and
 *     uncompressed size attributes. This shouldn't be a problem for x64
 *     processors.
 *
 *   - ZRAN_IMPORT_INCONSISTENT to indicate compressed size, or uncompressed
 *     size if known,   of the index file is inconsistent with the loaded
 *     compressed file.
 *
 *   - ZRAN_IMPORT_MEMORY_ERROR to indicate failure to allocate memory for new
 *     index. This typically result from out-of-memory.
 *
 *   - ZRAN_IMPORT_UNKNOWN_FORMAT to indicate given file is of unknown format.
 */
int zran_import_index(
  zran_index_t  *index, /* The index                  */
  FILE          *fd     /* Open handle to import file */
);

#endif /* __ZRAN_H__ */
