#!/bin/bash

echo "Start compiling..."

WORK_DIR="$(pwd)"
rm -rf "$WORK_DIR"/bin/*

mkdir -p "$WORK_DIR/bin/.o/"
cd "$WORK_DIR/bin/.o/"
for file in $(find "$WORK_DIR/src/c++/" -type f); do
    g++ -I "$WORK_DIR/include" -I "$WORK_DIR/lib" --std=c++11 -c "$file"
done

cd "$WORK_DIR"/bin/
g++ $(find $WORK_DIR/bin/.o/ -name "*.o" ! -name "main.o" -type f | xargs) "$WORK_DIR"/bin/.o/main.o -o main

exit 0