# MC-DCFP-comparison
Estimate R by A1: calculate the R function &amp; by A2: generate the cdf of failure times &amp; A3: multiple dependences

Main:
在原代码的基础上，更改了
(1)alpha,这里有个错误，之前取值0.01，改为1；
(2)NHPP,关于较大的T的选取，之前取值为1倍平均到达时间，经调试发现当倍数为4-5倍及以上时，结果与图4对应曲线一致性较好，现在取5倍平均到达时间。
Test:
在Main的基础上，更改了
(3)退化量x的计算式，之前的退化量是按照最新更新的退化率乘以当前时间tau计算的，目前的退化量是在更新时的退化量加上新退化率乘以delta_tau计算的。
