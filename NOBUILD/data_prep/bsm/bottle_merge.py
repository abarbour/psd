#! /usr/bin/python
#
# 29nov2006 --jrw  original

import sys, os, shutil

#sys.path.append("/usr/local/pbo/bin")
import bottle
#from pbo_exceptions import *

#########################################################################
class BottleMerge:

    #####################################################################
    def bottle_append(self, base_bottle, append_bottle, allow_gaps=True, bigendian=None):
        """
        Merge contents of base_bottle and append_bottle.  Typical use
        is that append_bottle is subsequent to base_bottle, but this is
        not enforced.  Combined output is written to base_bottle.
        """

        # special case - if base_bottle does not exist, then just
        # copy append_bottle to it.
        if not os.path.exists(base_bottle):
            shutil.copyfile(append_bottle, base_bottle)
        else:
            self.bottle_merge(base_bottle, append_bottle, base_bottle, allow_gaps, bigendian)

    #####################################################################
    def bottle_merge(self, b1name, b2name, outname, allow_gaps=True, bigendian=None):
        """
        Merge contents of b1name and b2name, writing results to outname.

        b1name must start at same time or before b2name.  If the files
        overlap, they will be merged only if the corresponding samples
        have MISSING values in one of the files.  If they do not overlap,
        b2name will be appended to b1name.

        If allow_gaps is True, missing data between the files will be
        provided using the MISSING value.

        The argument bigendian controls the formatting of individual
        data values.
        """

        b1 = bottle.Bottle(b1name)
        b1.read_header()
        b1.read_data()

        b2 = bottle.Bottle(b2name)
        b2.read_header()
        b2.read_data()
        
        if b1.file_metadata['interval'] != b2.file_metadata['interval']:
            raise PBO_FileContentsError("Can not merge bottle files with different sampling rates: %s=%g %s=%g" % (b1name, b1.file_metadata['interval'], b2name, b2.file_metadata['interval']))

        if b1.file_metadata['data_type'] != b2.file_metadata['data_type']:
            raise PBO_FileContentsError("Can not merge bottle files with different data types: %s=%s %s=%s" % (b1name, bottle.BTL_TYPE[b1.file_metadata['data_type']], b2name, bottle.BTL_TYPE[b2.file_metadata['data_type']]))

        if b1.file_metadata['id'] != b2.file_metadata['id']:
            raise PBO_FileContentsError("Can not merge bottle files with different IDs: %s=%d %s=%d" % (b1name, b1.file_metadata['id'], b2name, b2.file_metadata['id']))

        if b1.file_metadata['start'] > b2.file_metadata['start']:
            raise PBO_FileContentsError("%s does not start at or before %s" % (b1name, b2name))

        b1_end = b1.file_metadata['start'] + ((b1.file_metadata['num_pts'] - 1) * b1.file_metadata['interval'])
        outdata = []

        if b1_end <= b2.file_metadata['start']:
            # b1 completes, then b2 begins.
            for datum in b1.data:
                outdata.append(datum)
            if b1_end < b2.file_metadata['start']:
                if not allow_gaps:
                    raise PBO_FileContentsError("There is a gap between end of %s and beginning of %s; option set to not allow gaps" % (b1name, b2name))
                else:
                    num_missing = int((b2.file_metadata['start'] - b1_end - b1.file_metadata['interval']) / b1.file_metadata['interval'])
                    for i in range(num_missing):
                        outdata.append(bottle.BTL_MISSING[b1.file_metadata['data_type']])
            else:
                num_missing = 0
            for datum in b2.data:
                outdata.append(datum)
            out_num_pts = b1.file_metadata['num_pts'] + num_missing + b2.file_metadata['num_pts']

        else:
            # b1 starts first, but they overlap - see if b2 fills in gaps
            for datum in b1.data:
                outdata.append(datum)
            idx = int((b2.file_metadata['start'] - b1.file_metadata['start']) / b1.file_metadata['interval'])
            for datum in b2.data:
                if datum == bottle.BTL_MISSING[b1.file_metadata['data_type']]:
                    idx = idx + 1
                elif outdata[idx] == bottle.BTL_MISSING[b1.file_metadata['data_type']]:
                    outdata[idx] = datum
                    idx = idx + 1
                elif outdata[idx] == datum:
                    idx = idx + 1
                else:
                    raise PBO_FileContentsError("Trying to replace non-null data in %s with data from %s; outdata[%d]=%s datum=%s " % (b1name, b2name, idx, ("%" + bottle.BTL_TYPE_TO_PRINT[b1.file_metadata['data_type']]) % outdata[idx], ("%" + bottle.BTL_TYPE_TO_PRINT[b1.file_metadata['data_type']]) % datum))
            out_num_pts = b1.file_metadata['num_pts']

        out_start = b1.file_metadata['start']
        out_interval = b1.file_metadata['interval']
        out_data_type = b1.file_metadata['data_type']
        out_id = b1.file_metadata['id']

        # we've gotten everything we need from the two source bottle
        # files.  close them so that if we are overwriting one, it is safe
        del b1
        del b2

        if bigendian != None:
            bottle.write_bottle(outname, out_start, out_interval, out_num_pts, \
                    out_data_type, out_id, outdata, bigendian=bigendian)
        else:
            bottle.write_bottle(outname, out_start, out_interval, out_num_pts, \
                    out_data_type, out_id, outdata)


#########################################################################
# do something if invoked directly

if __name__ == '__main__':

    usage = "usage: %prog [-h] [-g] [-b] bottle1 bottle2 outfile\n" + \
            "\n" + \
            "   If not overlapping, will append bottle2 to end of\n" + \
            "   bottle1, creating new bottle file 'outfile'.\n" +\
            "\n" + \
            "   If overlapping, will examine bottle2 to see if it fills\n" + \
            "   in gaps in bottle1, and create new bottle file 'outfile'."
    from optparse import OptionParser
    parser = OptionParser(usage)
    parser.set_defaults(allow_gaps = True)
    parser.set_defaults(bigendian = False)
    parser.add_option("-g", "--gap-abort", \
            help="The '-g' switch will cause operation to abort " + \
            "if there is a gap between files.", \
            action="store_false", dest="allow_gaps")
    parser.add_option("-b", "--big-endian", \
            help="The '-b' switch will result in the output bottle " + \
            "file being written in big-endian " + \
            "(SPARC/Solaris, PowerPC/Mac) format rather than " + \
            "little-endian (Intel/Linux&Windows) format.", \
            action="store_true", dest="bigendian")
    (opts, args) = parser.parse_args()

    try:
        if len(args) != 3:
            parser.error("incorrect number of arguments")
        else:
            m = BottleMerge()
            m.bottle_merge(args[0], args[1], args[2], \
                    allow_gaps=opts.allow_gaps, bigendian=opts.bigendian)
    except:
        print "ERROR:", sys.exc_info()[1]
        raise

