#!/usr/bin/env bash
### Test1
echo "Running test1..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test.mp.txt \
    -E 0.01 -e 0.01 \
    > t/test_out/test1.out.txt \
    2>/dev/null
diff t/correct_out/test1.correct.out.txt t/test_out/test1.out.txt

