w5

[matrixx1,matrixy1] = CurveshorteningflowWITHFORCINGw5(numberofpointsinrho,numberofpointsintime,tau);
[matrixx2,matrixy2] = CurveshorteningflowWITHFORCINGw5(numberofpointsinrho*2,numberofpointsintime*2,tau);
[matrixx3,matrixy3] = CurveshorteningflowWITHFORCINGw5(numberofpointsinrho*4,numberofpointsintime*4,tau);
[matrixx4,matrixy4] = CurveshorteningflowWITHFORCINGw5(numberofpointsinrho*8,numberofpointsintime*8,tau);

radattau = sqrt(1+ (2*tau));

rad1 = sqrt(matrixx1(numberofpointsintime+1,2)^2 + (matrixy1(numberofpointsintime+1,2))^2);
rad2 = sqrt((matrixx2((numberofpointsintime*2)+1,2)^2 + (matrixy2((numberofpointsintime*2)+1,2))^2));
rad3 = sqrt((matrixx3((numberofpointsintime*4)+1,2)^2 + (matrixy3((numberofpointsintime*4)+1,2))^2));
rad4 = sqrt((matrixx4((numberofpointsintime*8)+1,2)^2 + (matrixy4((numberofpointsintime*8)+1,2))^2));


display("The radius at t= " + tau + "Should be = " + radattau);
display("For number of points in rho = " + numberofpointsinrho + " And for number of points in time = " + numberofpointsintime + ", The calculated radius = " + rad1 )
display("For number of points in rho * 2 = " + numberofpointsinrho*2 + " And for number of points in time * 2 = " + numberofpointsintime*2 + ", The calculated radius = " + rad2 )
display("For number of points in rho * 4 = " + numberofpointsinrho*4 + " And for number of points in time * 4 = " + numberofpointsintime*4 + ", The calculated radius = " + rad3)
display("For number of points in rho * 8 = " + numberofpointsinrho*8 + " And for number of points in time * 8 = " + numberofpointsintime*8 + ", The calculated radius = " + rad4)