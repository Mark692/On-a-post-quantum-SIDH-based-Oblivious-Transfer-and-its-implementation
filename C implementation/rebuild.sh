rm seedot
gcc -o seedot seedot.c lib/libsike.a
x-terminal-emulator -e "bash -c './seedot; bash;'"
