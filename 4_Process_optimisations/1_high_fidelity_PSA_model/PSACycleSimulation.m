function [ objectives, constraints ] = PSACycleSimulation( x, IsothermParams, material_property,x0, type, N)

% Retrieve process variables
% x(1) = P_0         Adsorption pressure [Pa]
% x(2) = ndot_0      Inlet molar flux [mol/s/m^2]
% x(3) = t_ads       Time of adsorption step [s]
% x(4) = alpha       Light product reflux ratio [-]
% x(5) = epsilon     Void fraction of the bed [-]
% x(6) = r_p         Radius of the pellets [m]
% x(7) = e_p         pellet porosity [-]

process_variables = [x(1), x(1)*x(2)/8.314/313.15,  x(3),  x(4),  x(5), x(6),x(7)] ;

try
    [objectives, constraints] = PSACycle(process_variables, IsothermParams, material_property,x0, type, N) ;
catch
    %warning('Problem using function.  Assigning a value of 0 for objectives and constraints vialations');
    objectives(1)  = 1e5 ;
    objectives(2)  = 1e5 ;
    constraints(1) = 1 ;
    constraints(2) = 1 ;
    constraints(3) = 1 ;
end

end

