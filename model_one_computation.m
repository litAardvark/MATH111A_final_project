function result = model_one_computation(w, d, t)
b = 6;
c = 9;
ln_offset = 0.6;
t_p = 6;
d_f = 10; %dampening factor

alpha = 8*cos(d*pi)+52;
week_factor = 1;
if(w<=8)
 week_factor = (1+(w/8));
else
   week_factor = (1+(w/16));
end
alpha = alpha*week_factor;
In = (alpha/b)*log(c*(t-ln_offset));
Out = exp((alpha/b)*((t-t_p)/d_f));
result = In - Out;
end