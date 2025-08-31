#!/bin/bash

DIR="data/*/"
FILES=$(find data/*/ -type f)
INTERVAL=1
ALGORITHMS="conn-I-pla cov-pla optimal-pla swing-filter slide-filter semi-optimal-pla semi-mixed-pla"

declare -A ERROR_MAPPER
ERROR_MAPPER["humi"]="0.92 2.3 4.6 6.9 9.2"
ERROR_MAPPER["temp"]="1.123 2.807 5.615 8.42 11.23"
ERROR_MAPPER["wind"]="0.32 0.8 1.6 2.4 3.2"
ERROR_MAPPER["EV-power"]="1.3 3.253 6.5 9.76 13.01"
ERROR_MAPPER["ping"]="21.05 52.63 105.27 157.905 210.54"
ERROR_MAPPER["mem-avail"]="0.3912 0.983 1.966 2.95 3.932"
ERROR_MAPPER["home-power"]="0.685 1.712 3.423 5.135 6.847"
ERROR_MAPPER["JPXJPY"]="170.25 425.625 851.25 1276.875 1702.5"
ERROR_MAPPER["EURUSD"]="0.0057 0.014 0.028 0.043 0.057"
ERROR_MAPPER["BTCUSD"]="196.64 491.6 983.21 1474.82 1966.42"
ERROR_MAPPER["ECG"]="0.0325 0.081 0.1627 0.244 0.3254"
ERROR_MAPPER["wafer"]="0.152 0.379 0.759 1.138 1.52"
ERROR_MAPPER["fungi"]="0.871 2.177 4.355 6.533 8.71"
ERROR_MAPPER["starlight-curves"]="0.081 0.2 0.407 0.61 0.814"
ERROR_MAPPER["bus-lat"]="0.002 0.005 0.011 0.017 0.022"
ERROR_MAPPER["bus-lng"]="0.0032 0.008 0.0163 0.024 0.032"

for ALGO in $ALGORITHMS; do
    for INDEX in {0..4}; do
        for FILE in $FILES; do
            if [[ $FILE == *"synthesis"* ]]; then
                continue
            fi

            BASE=$(basename $FILE)
            NAME="${BASE%.*}"
            read -a ERROR <<< ${ERROR_MAPPER["$NAME"]}

            sed "s|<DATA>|$FILE|g; s|<NAME>|$NAME|g; s|<INTERVAL>|$INTERVAL|g; s|<ERROR>|${ERROR[$INDEX]}|g" "conf/template/$ALGO.template" > "conf/$ALGO.json"
            for i in {0..9}; do
                echo ""
                echo "###############"
                echo "Compress $NAME with ERROR=${ERROR[$INDEX]}"
                echo "###############"
                echo ""
                script/run.sh conf/$ALGO.json

                if [[ $? -ne 0 ]]; then
                    echo "Compression Failed on index $INDEX!!!"
                    exit -1
                fi
            done
        done
    done
done

rm -rf out/compress/*
rm -rf out/decompress/*

exit 0