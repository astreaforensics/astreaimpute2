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

echo "Running test2..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test.mp.txt \
    -E 0.01 -e 0.01 -V \
    > t/test_out/test2.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test2.correct.out.txt t/test_out/test2.out.txt

echo "Running test3..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test.mp.txt \
    -E 0.01 -e 0.01 -V -B 3 -G 3 \
    > t/test_out/test3.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test3.correct.out.txt t/test_out/test3.out.txt

echo "Running test4..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test2.mp.txt \
    -E 0.01 -e 0.01 \
    > t/test_out/test4.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test4.correct.out.txt t/test_out/test4.out.txt

echo "Running test5..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test3.mp.txt \
    -E 0.01 -e 0.01 \
    > t/test_out/test5.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test5.correct.out.txt t/test_out/test5.out.txt

echo "Running test6..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test4.mp.txt \
    -E 0.01 -e 0.01 \
    > t/test_out/test6.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test6.correct.out.txt t/test_out/test6.out.txt

echo "Running test7..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test5.mp.txt \
    -E 0.01 -e 0.01 \
    > t/test_out/test7.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test6.correct.out.txt t/test_out/test7.out.txt

echo "Running test8..."
./astrea-impute2.14 \
    -H t/test.hap.txt \
    -m t/test6.mp.txt \
    -E 0.01 -e 0.01 \
    > t/test_out/test8.out.txt \
    2>/dev/null
diff -I '^##filedate' t/correct_out/test6.correct.out.txt t/test_out/test8.out.txt
