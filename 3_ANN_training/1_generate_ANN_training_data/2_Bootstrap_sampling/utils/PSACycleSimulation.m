function [ objectives, constraints ] = PSACycleSimulation( x, IsothermParams,type, N, x0)

% wrapper function to pass inputs into PSACycle.m function which simulates the 
% PSA cycle

process_variables = [x(1),x(1)*x(2)/8.314/313.15,x(3),x(4),x(5),x(6),x(7),x(8),x(9)] ;
IsothermPar = [IsothermParams(round(x(10)),:),0];

	try
    [objectives, constraints] = PSACycle(process_variables, IsothermPar, type, N, x0) ;
	catch
    %warning('Problem using function.  Assigning a value of 0 for objectives and constraints vialations');
    objectives(1)  = 1e5 ;
	objectives(2)  = 1e5 ;
	constraints(1) = 1 ;
	constraints(2) = 1 ;
	constraints(3) = 1 ;
	end
	
end

