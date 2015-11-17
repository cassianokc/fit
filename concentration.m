function c = concentration(t, fixedParams, params)
	c0 = fixedParams(1);
	res = params(1);
	k = params(2);
    c = res+(c0-res).*exp(-k.*t);
end
