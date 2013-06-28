Installing on Debian / Ubuntu
=============================

    $ sudo apt-get install fort77 xutils-dev libc-dev pvm-dev bison \
        libnetcdf-dev
    $ find . \( -name '*.o' -o -name '*.a' \) -exec rm {} \;
    $ rm -f Makefile mctwostream/Makefile
    $ ln -nfs Makefile.linux Makefile
    $ ln -nfs Makefile.linux mctwostream/Makefile
    $ (cd mctwostream; make depend && make libtwostream.a OBJDIR=.)
    $ make depend
    $ make optp
