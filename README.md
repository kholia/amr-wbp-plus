#### Basic Information

This code was downloaded from [here](http://www.3gpp.org/ftp/Specs/archive/26_series/26.304/)
and is based on `26304-f00.zip` release.

Floating-point implementation -> http://www.3gpp.org/ftp/Specs/archive/26_series/26.304/

Fixed-point implementation -> http://www.3gpp.org/ftp/Specs/archive/26_series/26.273/

3GPP file format -> http://www.3gpp.org/ftp/Specs/archive/26_series/26.412

I would love to port this to FFmpeg (assuming no licensing issues...) - looking
for FFmpeg mentors!

#### Usage

```
$ wget http://www.music.helsinki.fi/tmt/opetus/uusmedia/esim/a2002011001-e02-16kHz.wav

$ wine encoder.exe -rate 16 -if a2002011001-e02-16kHz.wav -of out.3gp
3gp File Format
Encoding @  16.00kbps
 --- Running ---

$ mediainfo out.3gp
General
Complete name                            : out.3gp
Format                                   : MPEG-4
Format profile                           : 3GPP Media Release 6 Basic
Codec ID                                 : 3gp6 (3gp6/isom)
File size                                : 110 KiB
Duration                                 : 54 s 160 ms
Overall bit rate mode                    : Constant
Overall bit rate                         : 16.7 kb/s
Encoded date                             : UTC 2017-02-10 06:10:58
Tagged date                              : UTC 2017-02-10 06:10:58

Audio
ID                                       : 1
Format                                   : sawp
Codec ID                                 : sawp
Duration                                 : 54 s 160 ms
Bit rate mode                            : Constant
Bit rate                                 : 16.2 kb/s
Channel(s)                               : 2 channels
Sampling rate                            : 6 464 Hz
Bit depth                                : 16 bits
Stream size                              : 107 KiB (97%)
Encoded date                             : UTC 2017-02-10 06:10:58
Tagged date                              : UTC 2017-02-10 06:10:58

$ wine decoder.exe -if out.3gp -of out.wav
...
File type brand '3gp6' minor version 256
Compatible brands: 'isom' (yes) '3gp4' (no) '3gp5' (no)
Track 1 'sawp'  Duration 54.160s  Time scale 72000
Decoding @  16.00kbps

$ mediainfo out.wav
General
Complete name                            : out.wav
Format                                   : Wave
File size                                : 9.92 MiB
Duration                                 : 54 s 160 ms
Overall bit rate mode                    : Constant
Overall bit rate                         : 1 536 kb/s

Audio
Format                                   : PCM
Format settings, Endianness              : Little
Format settings, Sign                    : Signed
Codec ID                                 : 1
Duration                                 : 54 s 160 ms
Bit rate mode                            : Constant
Bit rate                                 : 1 536 kb/s
Channel(s)                               : 2 channels
Sampling rate                            : 48.0 kHz
Bit depth                                : 16 bits
Stream size                              : 9.92 MiB (100%)
```

##### Links

- https://github.com/MPEGGroup/isobmff/
