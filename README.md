# cpp_language_feature

## Purpose
This project aims to check the cpp language features. 
(only `__builtin_assume` currently)

## How to run
### Get IR file
```bash
clang if_test.c -S -emit-llvm -O2 -o if_test.ll
```

### Get EXE file
```bash
clang if_test.c -O2 -o if_test
```
