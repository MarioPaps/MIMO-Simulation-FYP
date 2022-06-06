function y = findzero(b,c,x0)

y = fzero(@poly,x0);

   function y = poly(x)
   y = x^3 + b*x + c;
   end
end