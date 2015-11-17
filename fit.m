function [minElem, minErr] = fit(x, y, lowerBound, upperBound, fixedParams, func)
popSize = 100;
dimSize = length(lowerBound);

minErr = realmax;
minElem = zeros(1, dimSize);
err = zeros(1, popSize);
mutationFactor=0.3;
lambda=0.7;
recProb=0.8;


elems = (rand(popSize, dimSize)-lowerBound./upperBound).*upperBound;

for i=1:popSize
		calcErr = abs(y-func(x, fixedParams,elems(i,:)));
		err(1,i) = mean(calcErr.^2);
end

[minErr, minIdx] = min(err(1,:))
minElem = elems(minIdx,:)



for it = 1:100000;
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
	[minErr, minIdx] = min(err(1,:))
	minElem = elems(minIdx,:)
end

