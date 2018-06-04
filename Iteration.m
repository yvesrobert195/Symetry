addpath(genpath('../'))
%% User Input
fprintf('*********************************************************************\n')
fprintf('ORIFICE ITERATION TIME SCRIPT\n')
fprintf(datestr(now))
fprintf('\n*********************************************************************\n')
User_Input
% INPUT SPECIFICALLY FOR ITERATION TIME
comment=['_A']; %additional comments for the name (put '_')
sol_name=['sym_' datestr(now,'mm-dd-yy_HH-MM-SS') comment '.mat']; %creates unique mat file for the test
ini=10;
step=10;

%% Initialization
fprintf('*********************************************************************\n')
fprintf('Initialization\n')
fprintf('*********************************************************************\n')
% Assign geometry
fprintf('\tReading geometry\n')
Geometry=geometry(Input.Core);

% Reads det files
fprintf('\tReading power data\n')
[Pb.Constraints.same_pos,Input.Q,Input.adjacentAssemblies]= readQ_sym(Input.powerDetectorFiles,Input.assemblyPowerThreshold);
Input.Q_ave = sum(Input.Q,2)/length(Input.powerDetectorFiles); % divide by number of steps to get average assembly power over cycle

fprintf('\tManipulating data\n');
Pb.Var.nass = length(Input.Q_ave); % number of assemblies in problem
Pb.Var.nsteps = length(Input.powerDetectorFiles);



n=0;
while 1
    fprintf('\n\n*********************************************************************\n')
    fprintf(['TEST WITH ' num2str(ini+n*step) ' FLOWRATES\n'])
    fprintf('*********************************************************************\n')
    
    Pb.Var.x=logspace(-2,2,ini+n*step);
    Pb.Var.npossflows = length(Pb.Var.x); % number of possible flowrates specified as data
    Pb.Var.nvars = Pb.Var.nass+Pb.Var.nass*Pb.Var.npossflows+Pb.Var.npossflows; % number of total variables
    
    % Build tables of coolant properties and temperature
    fprintf('\tBuilding coolant properties and temperature tables\n')
    Coolant=coolant_properties('sodium',Input,Pb,Geometry);
    
    % Find which assemblies have adjacents and how many adjacent pairs
    Input.nadj = nnz(Input.adjacentAssemblies); % number of adjacent assembly pairs
    
    % Read and create rings
    Geometry.rings=find_rings(Input.adjacentAssemblies);
    
    % Create constraints
    fprintf('*********************************************************************\n')
    fprintf('Constraint matrices\n')
    fprintf('*********************************************************************\n')
    Pb=make_constraints(Pb,Input,Coolant,Geometry);
    
    % Solve
    fprintf('*********************************************************************\n')
    fprintf('Solving with CPLEX\n');
    fprintf('*********************************************************************\n')
    [Solution.solutionvector, Solution.objval, Solution.status, Solution.output] = cplexmilp(Pb.CPLEX.c, Pb.CPLEX.Aineq, Pb.CPLEX.bineq, Pb.CPLEX.Aeq, Pb.CPLEX.beq, [], [], [], [], [], Pb.CPLEX.ctype,[],Pb.CPLEX.opts);
    fprintf('exit status = % i\n', Solution.status);
    fprintf('solution time = % f\n', Solution.output.time);
    
    % Post processing and check
    if Solution.status==1
        fprintf('*********************************************************************\n')
        fprintf('SOLUTION FOUND : Post-processing and checking errors\n');
        fprintf('*********************************************************************\n')
        Output=post_process(Solution,Pb,Input,Coolant,Geometry);
    else
        Output={};
        fprintf('\n**********************    NO SOLUTION FOUND    **********************\n\n')
    end
    
    Results(n+1,:)={length(Pb.Var.x) Output Solution Solution.status};
    n=n+1;
    save(['../Solutions/' sol_name],'Results')
end


fprintf('*********************************************************************\n')
fprintf('END OF SCRIPT\n');
fprintf('*********************************************************************\n')
diary off
try
    movefile clone* ./Clones/;
end