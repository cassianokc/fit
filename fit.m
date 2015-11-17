function [minElem, minErr] = fit(x, y, lowerBound, upperBound, fixedParams, func)
%Finds the parameters of a function that minimizes the MSE of the obtained func results for a vector x with the expected/measured values y by using diferential evolution.
%	x: inputs of f where the MSE are being calculated
%	y: expected/measured outputs of the functions
%	lowerBound: vector containing the lowest possible values for all the
%parameters in the same order as they are used on func
%	upperBound: vector containing the highest possible values for all the
%parameters in the same order as they are used on func
%	fixedParams: these variable is always passed to func as a parameter without any changes, this allows more flexibility for defining the function
% function: A function handler/reference that is going to be called for
%computing the its values on input vector x with different parameters. The
%function must receive 3 parameters, the input vector x, the constant values
%fixedParams and a vector containing the parameters that are going to be fit
% return minElem: the parameters that gave the provide the minimum MSE
% return minErr: the minimum MSE found

popSize = 100;
dimSize = length(lowerBound);

minErr = realmax;
minElem = zeros(1, dimSize);
err = zeros(1, popSize);
mutationFactor=0.3;
recProb=0.8;


elems = (rand(popSize, dimSize)-lowerBound./upperBound).*upperBound;

for i=1:popSize
		calcErr = abs(y-func(x, fixedParams,elems(i,:)));
		err(1,i) = mean(calcErr.^2);
end

[minErr, minIdx] = min(err(1,:));
minElem = elems(minIdx,:);
it = 0;




while it<100
	for i=1:popSize
		%randomize the parents
		parentIdx = round(rand(1, 3)*(popSize-1))+1;
		while (unique([parentIdx(:); i]) != 4)
			parentIdx = round(rand(1, 3)*(popSize-1))+1;
		end
		% calculate the donor
		u = elems(parentIdx(1),:) + mutationFactor*(elems(parentIdx(2),:)-elems(parentIdx(3),:));
		u = elems(parentIdx(1),:).*(u<lowerBound || u>upperBound) + u.*(u>lowerBound && u<upperBound);
		r = rand(1, dimSize);
		intR = rand(1)*(dimSize-1)+1;
		%Get the new vector
		v = elems(i,:).*(r>=recProb && i != intR) + u.*(r<recProb |  i==intR);

		calcErr = abs(y-func(x, fixedParams, v));
		vErr = mean(calcErr.^2);
		%If it was better then previous value, update
		if vErr < err(1, i)
			err(1, i) = vErr;
			elems(i,:) = v;
		end
	end
	lastMinElem = minElem;
	[minErr, minIdx] = min(err(1,:));
	minElem = elems(minIdx,:);
	if max(max(lastMinIdx-minIdx)) < 0.0001
		it = it+1;
	end
end

