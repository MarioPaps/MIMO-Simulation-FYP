Ki=10;
M=2;
mod1= deg2rad(43);
mod2= deg2rad(2);
bitstream1= round(rand(1,Ki*2)); %random bitstream->i could use this in my function
bitstream2= round(rand(1,Ki*2));

m1= QPSKMod(bitstream1,sqrt(2),mod1);
m2= QPSKMod(bitstream2,0.5,mod2);