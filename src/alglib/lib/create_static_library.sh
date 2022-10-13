#!/bin/bash
wget http://www.alglib.net/translator/re/alglib-3.9.0.cpp.gpl.tgz
mv ./alglib-3.9.0.cpp.gpl.tgz /usr/local/alglib-3.9.0.cpp.gpl.tgz
cd /usr/local/
tar -xvzf alglib-3.9.0.cpp.gpl.tgz
mv ./cpp ./alglib
cd ./alglib
mv ./src ./include
mkdir ./lib
cd ./include
g++ -c *.cpp
ar rcs alglib.a *.o
rm -rf *.o
mv ./alglib.a ../lib/libalglib.a
echo "alglib Installed successfully!"
cd ~
