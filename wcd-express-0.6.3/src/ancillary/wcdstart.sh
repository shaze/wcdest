#! /bin/sh
### BEGIN INIT INFO
# Provides:          wcdstart.sh
# Required-Start:
# Required-Stop:
# Default-Start:     S
# Default-Stop:
# Short-Description:Sets up wcd.
# Description:
### END INIT INFO

PATH=/sbin:/bin:/usr/bin


do_start () {
    /bin/rm -rf /tmp/wcd*
    wget http://wcd-source.s3.amazonaws.com/wcd-express.tar.gz  -P /tmp
    wget http://wcd-source.s3.amazonaws.com/setup.sh -P /tmp
    /bin/rm -f /usr/local/bin/wcd*
    /bin/rm -f /usr/local/share/info/wcd*
    cd /tmp
    tar -xzf wcd-express.tar.gz
    cd wcd-express
    ./configure
    make
    make install
}

case "$1" in
  start|"")
        do_start
        ;;
  restart|reload|force-reload)
        echo "Error: argument '$1' not supported" >&2
        exit 3
        ;;
  stop)
        # No-op
        ;;
  *)
        echo "Usage: wcdstart.sh [start|stop]" >&2
        exit 3
        ;;
esac

