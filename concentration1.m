function c = concentration1(t, fixedParams, params)
	initialC1 = fixedParams(1);
	initialC2 = fixedParams(1);
	lastT = t(1);
	res1 = params(1);
	k1 = params(2);
	th = params(3);
	res2 = params(4);
	k2 = params(5);
	for i=1:length(t)
		if initialC2 > th
    		c(i) = res1+(initialC1-res1).*exp(-k1*t(i));
			initialC2 = c(i);
			lastT = t(i);
		else
    		c(i) = res2+(initialC2-res2).*exp(-k2*(t(i)-lastT));
		end
	end
	c = c';
end
