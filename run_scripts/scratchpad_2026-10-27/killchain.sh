pkill -f chainAB.sh; pkill -f 'ephproto.bandcheck'; pkill -f wvar.py; sleep 2; pgrep -fa 'chainAB|bandcheck|wvar' || echo clean
