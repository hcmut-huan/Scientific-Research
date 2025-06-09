#!/bin/bash

# Preprocessing config
python3 src/python/validate.py $1 > .temp

if [[ $? == 1 ]]; then
    echo "Algorithm is not supported."
    exit 1
elif [[ $? == 2 ]]; then
    echo "Missing required parameter."
    exit 2
elif [[ $? == 3 ]]; then
    echo "Invalid parameter."
    exit 3
fi

DATA=$(sed -n "1p" .temp)
COMPRESS=$(sed -n "2p" .temp)
DECOMPRESS=$(sed -n "3p" .temp)
INTERVAL=$(sed -n "4p" .temp)
ALGO=$(sed -n "5p" .temp)

mkdir -p out/compress
mkdir -p out/decompress
rm -f .temp

# Compressing phase
echo -e "-------------------------"
echo "Start compressing..."
echo '' > $COMPRESS
echo '' > $COMPRESS.mon
echo '' > $COMPRESS.time
bin/compress $DATA $COMPRESS $ALGO

# Decompressing phase
echo -e "\n-------------------------"
echo "Start decompressing..."
echo '' > $DECOMPRESS
echo '' > $DECOMPRESS.mon
echo '' > $DECOMPRESS.time
bin/decompress $COMPRESS $DECOMPRESS $INTERVAL $ALGO

# Statistic phase
echo -e "\n-------------------------"
echo "Start statisticizing..."
python3 src/python/statistics.py $DATA $DECOMPRESS $COMPRESS > .statistic

echo -n $DATA,$(echo $ALGO | awk -F " " '{print $1}'),$(echo $ALGO | awk -F " " '{print $2}') >> out/experiments.csv
cat .statistic | while read line; do
    echo $line
    echo -n ,$(echo $line | awk -F ":" '{print $2}' | xargs) >> out/experiments.csv 
done

echo -n ,$(cat $COMPRESS.time | grep -oE '[0-9]+\.[0-9]+|[0-9]+' | paste -sd, -) >> out/experiments.csv
echo ,$(cat $DECOMPRESS.time | grep -oE '[0-9]+\.[0-9]+|[0-9]+' | paste -sd, -) >> out/experiments.csv
rm -f .statistic

exit 0
