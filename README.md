# Grasshopper
This project is optimized implementation of Russian ГОСТ Р 34.12-2015 Kuznechik (Grasshopper).  
## Results:  
ECB encryption, Mb/s

|                 | i3-4030u | i5-5200u | i5-7200u | i5-8250u |
|-----------------|----------|----------|----------|----------|
| Baseline        | 6        | 8.8      | 10.2     | 11.4     |
| With LUT        | 53       | 77       | 99       | 107      |
| LUT, vector xor | 69       | 99       | 131      | 139      |
| Offset precalc  | 79       | 113      | 135      | 142      |
| AVX2            | 84       | 121      | 150      | 159      |

CFB decryption, Mb/s

|                 | i3-4030u | i5-5200u | i5-7200u | i5-8250u |
|-----------------|----------|----------|----------|----------|
| Baseline        | 5.2      | 8.7      | 10.2     | 9.4      |
| With LUT        | 50       | 73       | 94       | 103      |
| LUT, vector xor | 73       | 105      | 131      | 145      |
| Offset precalc  | 80       | 118      | 135      | 150      |
| AVX2            | 82       | 120      | 145      | 159      |
