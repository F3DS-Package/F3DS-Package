#!/bin/bash
PREFIX=$1

cat <<EOF > $PREFIX/setenv.sh
F3DS_ROOT=$PREFIX; export F3DS_ROOT
EOF
cat <<'EOF' >> $PREFIX/setenv.sh
PATH=$PATH:$F3DS_ROOT/bins; export PATH
F3DS_LIBS=$F3DS_ROOT/libs; export F3DS_LIBS
F3DS_MODS=$F3DS_ROOT/mods; export F3DS_MODS
EOF