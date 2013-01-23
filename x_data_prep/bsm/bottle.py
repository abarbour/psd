#! /usr/bin/python        
#
# 12jul2006 --jrw  original

import sys, os, datetime, struct

#sys.path.append("/usr/local/pbo/bin")
#from pbo_exceptions import *

BTL_EPOCH = datetime.datetime(1970, 1, 1, 0, 0, 0, 0, None)

BTL_TYPE = {0: 'INT2', 1: 'INT4', 2: 'FLOAT4'}
BTL_MISSING = {0: 32767, 1: 999999, 2: 999999}
BTL_TYPE_SIZE = {0: 2, 1: 4, 2: 4}
BTL_TYPE_TO_PACK = {0: 'h', 1: 'i', 2: 'f' }
BTL_TYPE_TO_PRINT = {0: 'd', 1: 'd', 2: 'g' }

#########################################################################
def write_bottle(filename, start, interval, num_pts, data_type, id, data, bigendian=False):
    """
    """

    if num_pts != len(data):
        raise PBO_FileContentsError("When writing bottle file %s - num_pts=%d, len(data)=%d" % (filename, num_pts, len(data)), filename)

    if data_type not in [0, 1, 2]:
        raise PBO_FileContentsError("When writing bottle file %s - invalid data type %d" % (filename, data_type), filename)

    if bigendian:
        endian = '>'
    else:
        endian = '<'
    magic = 0x015D
    unused = 0
    header_size = 40
    usgs_lock = 0
    header = struct.pack("%shhidfiiiii" % endian, magic, unused, \
                    header_size, start, interval, num_pts, \
                    data_type, BTL_MISSING[data_type], usgs_lock, id)
    if len(header) != header_size:
        raise PBO_FileContentsError("Invalid header size?!?!", filename)

    file = open(filename, 'wb')
    file.write(header)
    for datum in data:
        bytes = struct.pack(endian + BTL_TYPE_TO_PACK[data_type], datum)
        file.write(bytes)
    file.close()

#########################################################################
class Bottle:
    def __init__(self, filename):
        """
        """

        self.file_metadata = {}
        self.file_metadata['filename'] = filename
        if not os.path.isfile(self.file_metadata['filename']):
            raise IOError, "Can not find file '%s'" % self.file_metadata['filename']

        self.file = open(self.file_metadata['filename'], 'rb')
        magic = self.file.read(2)
        if magic[0] == chr(0x01) and magic[1] == chr(0x5D):
            self.file_metadata['endian'] = '>'
        elif magic[0] == chr(0x5D) and magic[1] == chr(0x01):
            self.file_metadata['endian'] = '<'
        else:
            print magic
            print "Not recognized as a bottle file: %s" % self.file_metadata['filename']
            raise PBO_FileContentsError("Not a bottle", self.file_metadata['filename'])

    ##################################################################
    def parse_filename(self):
        """
        Return info encoded within file name of bottle file from GTSM
        Logger.  Return is (stn4char, channel, year, doy, hour, min)
        where hour and min may be None.  All values are strings, even
        if they contain only digits.
        """

        # FIXME: we are hardcoding here knowledge of how to pick apart
        # a GTSM file name.  no other option really, but just be wary
        # of trusting file names.  and if you feed in something which
        # isn't a GTSM bottle file name, you'll get junk back - or
        # raise an exception.

        f = self.file_metadata['filename']

        stn4char = f[0:4]
        year = "20" + f[4:6]
        dayofyear = f[6:9]

        if f[-3:] == '_20' and len(f) == 19:
            # this is a Min file
            hour = f[9:11]
            min = f[11:13]
            channel = f[13:16]

        elif f[-3:] == '_20' and len(f) == 17:
            # this is the concatenation of 60 Min files
            hour = f[9:11]
            min = None
            channel = f[11:14]

        elif f[-3:-1] == 'CH' and len(f) == 14:
            # this is an Hour file
            hour = f[9:11]
            min = None
            channel = f[11:14]

        else:
            # well, we have to assume it is a Day file
            hour = None
            min = None
            channel = f[9:]

        return (stn4char, channel, year, dayofyear, hour, min)

    #####################################################################
    def read_header(self, print_it=False):
        """
        """

        self.file.seek(0)
        data = self.file.read(40)
        format = "%shhidfiiiii" % self.file_metadata['endian']
        (   self.file_metadata['magic'], \
            self.file_metadata['unused'], \
            self.file_metadata['header_size'], \
            self.file_metadata['start'], \
            self.file_metadata['interval'], \
            self.file_metadata['num_pts'], \
            self.file_metadata['data_type'], \
            self.file_metadata['missing'], \
            self.file_metadata['usgs_lock'], \
            self.file_metadata['id'] \
        ) = struct.unpack(format, data)
        self.file_metadata['start_timestamp'] = datetime.datetime.isoformat(BTL_EPOCH + datetime.timedelta(0, self.file_metadata['start']), ' ')

        if print_it:
            print ""
            print "file:        %s" % self.file_metadata['filename']
            print "magic:       %x" % self.file_metadata['magic']
            print "unused:      %d" % self.file_metadata['unused']
            print "header_size: %d" % self.file_metadata['header_size']
            print "start:       %s" % self.file_metadata['start_timestamp']
            print "interval:    %g" % self.file_metadata['interval']
            print "num_pts:     %d" % self.file_metadata['num_pts']
            if self.file_metadata['data_type'] in BTL_TYPE.keys():
                print "data_type:   %s" % BTL_TYPE[self.file_metadata['data_type']]
            else:
                print "data_type:   %d" % self.file_metadata['data_type']
            print "missing:     %d" % self.file_metadata['missing']
            print "usgs_lock:   %d" % self.file_metadata['usgs_lock']
            print "id:          %d" % self.file_metadata['id']
            print ""

        return self.file_metadata

    #####################################################################
    def read_data(self, print_it=False, timestamps=True):
        """
        """


        if print_it:
            if timestamps:
                print_format = "%%s %%%s" % BTL_TYPE_TO_PRINT[self.file_metadata['data_type']]
                # We use the python timedelta class, since it can
                # better represent intervals such as 0.05.  this becomes
                # (0,0,500000) rather than the floating point approximation
                # of 0.05 which prints as 0.050000000000000003.
                seconds = int(self.file_metadata['interval'])
                microseconds = int((self.file_metadata['interval'] - seconds) * 1000000)
                interval = datetime.timedelta(0, seconds, microseconds)
                start_offset = datetime.timedelta(0, self.file_metadata['start'])
            else:
                print_format = "%%%s" % BTL_TYPE_TO_PRINT[self.file_metadata['data_type']]

        pack_format = "%s%s" % (self.file_metadata['endian'], BTL_TYPE_TO_PACK[self.file_metadata['data_type']])
        self.file.seek(self.file_metadata['header_size'])
        self.data=[]
        for i in range(self.file_metadata['num_pts']):
            datum = self.file.read(BTL_TYPE_SIZE[self.file_metadata['data_type']])
            if datum == '':
                print "error: read() returned empty string?"
                continue
            if len(datum) != BTL_TYPE_SIZE[self.file_metadata['data_type']]:
                print "error: read() returned unexpected number of bytes"
                continue
            value = struct.unpack(pack_format, datum)
            self.data.append(value[0])
            if print_it:
                if timestamps:
                    print print_format % (datetime.datetime.isoformat(BTL_EPOCH + start_offset + (i * interval), ' '), value[0])
                else:
                    print print_format % (value[0])

        return self.data

#########################################################################
# do something if invoked directly

if __name__ == '__main__':
    def usage():
        print "Usage: %s [ -h | -d | -t ] filename" % sys.argv[0]
        print "      default: print header and all data with timestamps"
        print "       -h      print header info only"
        print "       -d      print data only, no header, no timestamps"
        print "       -t      print data only, no header, with timestamps"

    if len(sys.argv) == 2:
        b = Bottle(sys.argv[1])
        b.read_header(print_it=True)
        b.read_data(print_it=True)
    elif len(sys.argv) == 3 and sys.argv[1] == '-h':
        b = Bottle(sys.argv[2])
        b.read_header(print_it=True)
    elif len(sys.argv) == 3 and sys.argv[1] == '-d':
        b = Bottle(sys.argv[2])
        b.read_header(print_it=False)
        b.read_data(print_it=True, timestamps=False)
    elif len(sys.argv) == 3 and sys.argv[1] == '-t':
        b = Bottle(sys.argv[2])
        b.read_header(print_it=False)
        b.read_data(print_it=True)
    else:
        usage()

