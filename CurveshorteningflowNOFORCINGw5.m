function [matrixx,matrixy] = CurveshorteningflowNOFORCINGw5(numberofpointsinrho,numberofpointsintime,tau)

        % We look to implement Curve shortening flow using finite elements where we
        % have parametrised some intial curve as vectors of (x,y). This will be for
        % the unit circle.

                              %   Initial conditions
                                
        % For \rho in [0,1] we can parameterise the unit circle as
        % [cos2\pi\rho, sin2\pi\rho]^{T} depending on how many angle segments we choose
        % to implement. Since this is a closed unit circle, it has periodic
        % boundary conditions. We will include a forcing term (f) so that it
        % can be used such that forcing is included.


        % Create the mesh for rho in [0,1] and time in [0,tau] where tau is the end
        % point (note as tau -> 0.5 we enter problems for p = 0 (no
        % forcing))

        drho = (1-0)/numberofpointsinrho;
        dt = (tau-0)/numberofpointsintime;
 
        pointrho = zeros(1,numberofpointsinrho+1);
        pointtime = zeros(1,numberofpointsintime+1);
        
            for i = 1:numberofpointsinrho+1
                pointrho(i) = (i-1)*drho;
            end
        
            for j=1:numberofpointsintime+1
                pointtime(j) = (j-1)*dt;
            end

        % First create matrices for the x and y component of each vector

        matrixx = zeros(numberofpointsintime+1,numberofpointsinrho+1); % x component 
        matrixy = zeros(numberofpointsintime+1,numberofpointsinrho+1); % y component
        
        % Create column vectors for RHS linear coefficients of prev
        % iteration
        
        RHSCIvec = zeros(numberofpointsinrho,1);
        
        % Create Column vector for forcing term with coefficients
        
        RHSforcingvecx = zeros(numberofpointsinrho,1);
        RHSforcingvecy = zeros(numberofpointsinrho,1);

        % Set up initial conditions into the vector matrices above with
        % matrixx_0 = cos(2*pi*rho), matrixy_0 = sin(2*pi*rho).
        
        for i = 1:numberofpointsinrho+1
            matrixx(1,i) = cos(2*pi*pointrho(i));
            matrixy(1,i) = sin(2*pi*pointrho(i));
        end
        
        % We now look to build all of our matrix systems. Since these
        % change on each loop, we now start our iterated loop sequence over
        % time.
        
        j=1;
        
        
        
            for j=1:numberofpointsintime
                
            %-------------------------------------------------------------%
            
                    % Forcing term to start we will say
        
                    p = @(t) 0;
        
                    % If we want to use the condition p = 2/ r(t) where r(t) denotes
                    % the radius at time t , then we need to insert this into the loop so this
                    % will change every iteration.
                    
                    % How to calculate the radius, we will just use the first point (rho = 0) since the forcing is constant and
                    % we are starting with the unit circle where the curvature is constant.
                   
            
                    %Start with building the left hand side matrix to be
                    %inverted
                    
                    matrixtobeinverted = zeros(numberofpointsinrho,numberofpointsinrho);
                    
                    for i=2:numberofpointsinrho-1
                        
                    % Define the q vectors
                    q(i-1) = sqrt((matrixx(j,i) - matrixx(j,i-1))^2 + (matrixy(j,i) - matrixy(j,i-1))^2);
                    q(i) = sqrt((matrixx(j,i+1) - matrixx(j,i))^2 + (matrixy(j,i+1) - matrixy(j,i))^2);
            
                    % Define the coefficients of each term first
                    
                    % ------   LHS   ------
                    
                    LHcIminus1 = -1/(q(i-1));
                    LHCI =(1/(2*dt))*(q(i-1) + q(i)) + 1/(q(i-1)) + 1/(q(i));
                    LHCIplus1 = -1/(q(i));
            
                    % -----    RHS  ------
                    
                    RHCI(i) = (1/(2*dt))*(q(i-1) + q(i));
                    
                    % ------ RHS forcing term Coefficients term ----- %
                    
                    % First perpindicular term
                    
                    xnonperp1 = matrixx(j,i) - matrixx(j,i-1);
                    ynonperp1 = matrixy(j,i) - matrixy(j,i-1);
                    xnonperp2 = matrixx(j,i+1) - matrixx(j,i);
                    ynonperp2 = matrixy(j,i+1) - matrixy(j,i);
                    
                    % Using the perp operator in the clock wise direction
                    
                    perp1 = [ynonperp1;-xnonperp1];
                    perp2 = [ynonperp2;-xnonperp2];
                    
                    RHforcingterm = 0.5*(perp1 + perp2); % NOTE THIS IS A VECTOR
                    
                    matrixtobeinverted(i,i-1) = LHcIminus1 ;
                    matrixtobeinverted(i,i) = LHCI ;
                    matrixtobeinverted(i,i+1)= LHCIplus1 ;
                    
                    RHSCIvec(i) = RHCI(i) ;
                    
                    RHSforcingvecx(i) = dot(RHforcingterm,[1;0]);
                    RHSforcingvecy(i) = dot(RHforcingterm,[0;1]);
           
                    end
            
                    % Boundary conditions for LHS matrix --------------- %
                    
                    % using the periodicity property here
                    qb = sqrt((matrixx(j,1) - matrixx(j,numberofpointsinrho))^2 + (matrixy(j,1) - matrixy(j,numberofpointsinrho))^2);
               
                    matrixtobeinverted(1,numberofpointsinrho) = -1/(qb) ;
                    matrixtobeinverted(1,1) = (1/(2*dt))*(q(1) + (qb)) + 1/(q(1)) + 1/((qb)); 
                    matrixtobeinverted(1,2) = -1/(q(1));

                    % ---------------------------------------------- %
                    
                    matrixtobeinverted(numberofpointsinrho,numberofpointsinrho-1) = -1/q(numberofpointsinrho-1);
                    matrixtobeinverted(numberofpointsinrho,numberofpointsinrho) = (1/(2*dt))*(q(numberofpointsinrho-1) + qb) + 1/(q(numberofpointsinrho-1)) + 1/(qb);
                    matrixtobeinverted(numberofpointsinrho,1) = -1/(qb);
                    
                    % ---------------------------------------------- %
                    
                    % Boundary condition for RHS vector RHCI
                    
                    RHSCIvec(1) = (1/(2*dt))*(qb + q(1)); 
                    RHSCIvec(numberofpointsinrho) = (1/(2*dt))*(q(numberofpointsinrho-1) + qb);  

                    % Boundary condition for RHS forcing vector
                    
                    RHSforcingvecx(1) = dot(0.5*([(matrixy(j,1) - matrixy(j,numberofpointsinrho));-(matrixx(j,1) - matrixx(j,numberofpointsinrho))] + [(matrixy(j,2) - matrixy(j,1));-(matrixx(j,2) - matrixx(j,1))]),[1;0]);
                    RHSforcingvecx(numberofpointsinrho) = dot(0.5*([(matrixy(j,numberofpointsinrho) - matrixy(j,numberofpointsinrho-1));-(matrixx(j,numberofpointsinrho) - matrixx(j,numberofpointsinrho-1))] + [(matrixy(j,1) - matrixy(j,numberofpointsinrho));-(matrixx(j,1) - matrixx(j,numberofpointsinrho))]),[1;0]) ;
                    
                    RHSforcingvecy(1) = dot(0.5*([(matrixy(j,1) - matrixy(j,numberofpointsinrho));-(matrixx(j,1) - matrixx(j,numberofpointsinrho))] + [(matrixy(j,2) - matrixy(j,1));-(matrixx(j,2) - matrixx(j,1))]),[0;1]);
                    RHSforcingvecy(numberofpointsinrho) = dot(0.5*([(matrixy(j,numberofpointsinrho) - matrixy(j,numberofpointsinrho-1));-(matrixx(j,numberofpointsinrho) - matrixx(j,numberofpointsinrho-1))] + [(matrixy(j,1) - matrixy(j,numberofpointsinrho));-(matrixx(j,1) - matrixx(j,numberofpointsinrho))]),[0;1]) ;
                    
                    %-------------------------------------------------------------%
            
                    % First we solve for X - note only thing we need to be
                    % careful of is the forcing term, since the rest is
                    % constant
                    
                    invertedmatrix = inv(matrixtobeinverted);
                    newvecx = zeros(numberofpointsinrho,1);
                    
                    for i=1:numberofpointsinrho
                        newvecx(i) = matrixx(j,i) * RHSCIvec(i);
                    end
                    
                    newiterationx = invertedmatrix*(newvecx + p(pointtime(j))*RHSforcingvecx);
                    
                    for i = 1:numberofpointsinrho
                        matrixx(j+1,i) = newiterationx(i);
                    end
                    
                    % Remembering to set the periodic condition
                    
                    matrixx(j+1,numberofpointsinrho+1) = matrixx(j+1,1);
                    
                    %-----------------------------------------------------------%
                    % Second we solve for Y - note only thing we need to be
                    % careful of is the forcing term, since the rest is
                    % constant
                    newvecy = zeros(numberofpointsinrho,1);
                    
                    for i=1:numberofpointsinrho
                        newvecy(i) = matrixy(j,i) * RHSCIvec(i);
                    end

                    newiterationy = invertedmatrix*(newvecy + p(pointtime(j))*RHSforcingvecy);
                    
                    for i = 1:numberofpointsinrho
                        matrixy(j+1,i) = newiterationy(i);
                    end
            
                    % Remembering to set the periodic condition
                    
                    matrixy(j+1,numberofpointsinrho+1) = matrixy(j+1,1);
            
            end
        
        plot(matrixx(1,:),matrixy(1,:));
        hold on
        plot(matrixx(numberofpointsintime/2,:),matrixy(numberofpointsintime/2,:));
        hold on
        plot(matrixx(numberofpointsintime+1,:),matrixy(numberofpointsintime+1,:));

end
