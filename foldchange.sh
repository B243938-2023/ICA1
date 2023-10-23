#!/usr/bin/bash

awk '{
    if ($1 != "average") {
        if ($2 == 0) {
            foldchange = $1 / ($2 + 1e-10);
        } else {
            foldchange = $1 / $2;
        }
        log2fc = log(foldchange) / log(2);
        print log2fc;
    }
}'