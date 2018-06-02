function P=make_constraints_symetry(P,I,C,G)
nass=P.Var.nass;
npossflows=P.Var.npossflows;
nsteps=P.Var.nsteps;
nvars=P.Var.nvars;
x=P.Var.x;
Constraints.same_pos=P.Constraints.same_pos;
nadj=I.nadj;

Q=I.Q;

cp=C.heat_capacity;
rho=C.density;
Omega=C.T_gradient;
T_in=C.T_inlet;

% Constraint data
v_max = P.Constraints.v_max;
dT_max = P.Constraints.dT_max;
dP_max = P.Constraints.dP_max;
T_out_bar=P.Constraints.T_out_bar;
T_out_bar_tol=P.Constraints.T_out_bar_tol;

A_flow=G.Assembly.flow_area;

fprintf('\tInitialization of constraint matrices\n')
%find how many constraints
nineqs = nass*nsteps + 2*nsteps + nass*nsteps + 2*nadj*nsteps + npossflows*nass + nsteps*nass; %number of inequality constraints
neqs = nass + nass + 5*length(Constraints.same_pos); %number of equality constraints

%initialize constraint matrices and vectors
nelements_Aineq = nass*nsteps*npossflows+ 2*nsteps*nass*npossflows + nass*nsteps + 2*nsteps*nadj*npossflows + nass*npossflows + nsteps*nass*npossflows;
Aineq_i = zeros(nelements_Aineq,1);
Aineq_j = zeros(nelements_Aineq,1);
Aineq_v = zeros(nelements_Aineq,1);
bineq = zeros(nineqs, 1);
nelements_Aeq = (npossflows+1)*nass + npossflows*nass + 10*length(Constraints.same_pos);
Aeq_i = zeros(nelements_Aeq,1);
Aeq_j = zeros(nelements_Aeq,1);
Aeq_v = zeros(nelements_Aeq,1);
beq = zeros(neqs, 1);
P.CPLEX.c = zeros(nvars,1);

P.CPLEX.ctype = '';
for i = 1:nvars
    if i <= nass %mdot variables
        P.CPLEX.ctype(end+1) = 'C';
    elseif i <= nass+nass*npossflows %delta variables
        P.CPLEX.ctype(end+1) = 'B';
    else %beta variables
        P.CPLEX.ctype(end+1) = 'C';
    end
end

%%%%%
%specify the constraint matrices
%%%%%

fprintf('\tBuilding constraint matrices\n');

idx_Aineq = 1; %keep track of the number of elements in Aineq matrix
idx_Aeq = 1;

%% INEQUALITIES
%max outlet temp (constraint 1)
for k = 1:nsteps
    for i = 1:nass
        for j=1:npossflows
            Aineq_i(idx_Aineq) = (k-1)*nass+i;
            Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
            Aineq_v(idx_Aineq) = -x(j)*cp(i,j,k);
            idx_Aineq = idx_Aineq+1;
        end
        bineq((k-1)*nass+i) = -Q(i,k)/dT_max;
    end
end

%mixed outlet temp (constraint 2)
for k = 1:nsteps
    for i = 1:nass
        for j=1:npossflows
            Aineq_i(idx_Aineq) = nsteps*nass+k;
            Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
            if G.rings(i)<=P.Constraints.rings_outlet
                Aineq_v(idx_Aineq) = -x(j)*(Omega(i,j,k)+T_in-(T_out_bar-T_out_bar_tol));
            else
                Aineq_v(idx_Aineq) = 0;
            end
            idx_Aineq = idx_Aineq+1;
            
            Aineq_i(idx_Aineq) = nsteps*nass+nsteps+k;
            Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
            if G.rings(i)<=P.Constraints.rings_outlet
                Aineq_v(idx_Aineq) = x(j)*(Omega(i,j,k)+T_in-(T_out_bar+T_out_bar_tol));
            else
                Aineq_v(idx_Aineq) = 0;
            end
            idx_Aineq = idx_Aineq+1;
        end
    end
    bineq(nsteps*nass+k) = 0;
    bineq(nsteps*nass+nsteps+k) = 0;
end

%max flow (constraint 3)
for k=1:nsteps
    for i = 1:nass
        for j=1:npossflows
            Aineq_i(idx_Aineq) = nsteps*nass+2*nsteps+(k-1)*nass+i;
            Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
            Aineq_v(idx_Aineq) = x(j)-v_max*rho(i,j,k)*A_flow;
            idx_Aineq = idx_Aineq+1;
        end
    end
end
bineq(nsteps*nass+2*nsteps+1:nsteps*nass+2*nsteps+nsteps*nass) = 0;

%adjacent outlet temp (constraint 4)
constraint_idx = 1;
for k = 1:nsteps
    for i = 1:nass
        for ip = I.adjacentAssemblies(i,:)
            if ip > 0
                for j=1:npossflows
                    Aineq_i(idx_Aineq)=nsteps*nass+2*nsteps+nsteps*nass+constraint_idx;
                    Aineq_j(idx_Aineq)=nass+(i-1)*npossflows+j;
                    Aineq_v(idx_Aineq) = Omega(i,j,k);
                    idx_Aineq = idx_Aineq+1;
                end
                for j=1:npossflows
                    Aineq_i(idx_Aineq)= nsteps*nass+2*nsteps+nsteps*nass+constraint_idx;
                    Aineq_j(idx_Aineq)= nass+(ip-1)*npossflows+j;
                    Aineq_v(idx_Aineq)= -Omega(ip,j,k);
                    idx_Aineq = idx_Aineq+1;
                end
                if G.rings(i)<=G.nrings && G.rings(ip)<=G.nrings
                    bineq(nsteps*nass+2*nsteps+nsteps*nass+constraint_idx)=P.Constraints.xi_power;
                else
                    bineq(nsteps*nass+2*nsteps+nsteps*nass+constraint_idx)=P.Constraints.xi_blanket;
                end
                constraint_idx = constraint_idx + 1;
                
                for j=1:npossflows
                    Aineq_i(idx_Aineq)=nsteps*nass+2*nsteps+nsteps*nass+constraint_idx;
                    Aineq_j(idx_Aineq)=nass+(i-1)*npossflows+j;
                    Aineq_v(idx_Aineq) = -Omega(i,j,k);
                    idx_Aineq = idx_Aineq+1;
                end
                for j=1:npossflows
                    Aineq_i(idx_Aineq)= nsteps*nass+2*nsteps+nsteps*nass+constraint_idx;
                    Aineq_j(idx_Aineq)= nass+(ip-1)*npossflows+j;
                    Aineq_v(idx_Aineq)= Omega(ip,j,k);
                    idx_Aineq = idx_Aineq+1;
                end
                if G.rings(i)<=G.nrings && G.rings(ip)<=G.nrings
                    bineq(nsteps*nass+2*nsteps+nsteps*nass+constraint_idx)=P.Constraints.xi_power;
                else
                    bineq(nsteps*nass+2*nsteps+nsteps*nass+constraint_idx)=P.Constraints.xi_blanket;
                end
                constraint_idx = constraint_idx + 1;
            end
        end
    end
end

%number of groups (constraint 8)
for i = 1:nass
    for j = 1:npossflows
        Aineq_i(idx_Aineq) = nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+(i-1)*npossflows+j;
        Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
        Aineq_v(idx_Aineq) = 1;
        idx_Aineq = idx_Aineq+1;
        
        Aineq_i(idx_Aineq) = nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+(i-1)*npossflows+j;
        Aineq_j(idx_Aineq) = nass+nass*npossflows+j;
        Aineq_v(idx_Aineq) = -1;
        idx_Aineq = idx_Aineq+1;
    end
end
bineq(nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+1:nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+nass*npossflows) = 0;

% %Pressure loss (constraint 9)
for k=1:nsteps
    for i = 1:nass
        for j = 1:npossflows
            Aineq_i(idx_Aineq) = nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+nass*npossflows+(k-1)*nass+i;
            Aineq_j(idx_Aineq) = nass+(i-1)*npossflows+j;
            Aineq_v(idx_Aineq) = C.P_gradient(i,j,k);
            idx_Aineq = idx_Aineq+1;
        end
    end
end
bineq(nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+nass*npossflows+1:nsteps*nass+2*nsteps+nsteps*nass+2*nadj*nsteps+nass*npossflows+nass*nsteps) = P.Constraints.dP_max;

%% EQUALITIES
%select flow from discretized flows (constraint 5)
for i = 1:nass
    Aeq_i(idx_Aeq) = i;
    Aeq_j(idx_Aeq) = i;
    Aeq_v(idx_Aeq) = -1;
    idx_Aeq = idx_Aeq+1;
    for j = 1:npossflows
        Aeq_i(idx_Aeq) = i;
        Aeq_j(idx_Aeq) = nass+(i-1)*npossflows+j;
        Aeq_v(idx_Aeq) = x(j);
        idx_Aeq = idx_Aeq+1;
    end
    beq(i) = 0;
end

%select only one flowrate for a channel (constraint 6)
for i = 1:nass
    for j = 1:npossflows
        Aeq_i(idx_Aeq) = nass+i;
        Aeq_j(idx_Aeq) = nass+(i-1)*npossflows+j;
        Aeq_v(idx_Aeq) = 1;
        idx_Aeq = idx_Aeq+1;
    end
    beq(nass+i) = 1;
end

%Symetry
for i=1:length(Constraints.same_pos)
    for j=1:5
        Aeq_i(idx_Aeq) = nass+nass+5*(i-1)+j;
        Aeq_j(idx_Aeq) = Constraints.same_pos(1,i);
        Aeq_v(idx_Aeq) = 1;
        idx_Aeq = idx_Aeq+1;
        
        Aeq_i(idx_Aeq) = nass+nass+5*(i-1)+j;
        Aeq_j(idx_Aeq) = Constraints.same_pos(j+1,i);
        Aeq_v(idx_Aeq) = -1;
        idx_Aeq = idx_Aeq+1;
        
        beq(nass+nass+5*(i-1)+j) = 0;
    end
end

fprintf('\tBuilding the constraint matrices\n');

%build sparse matrices for Aineq and Aeq
P.CPLEX.Aineq = sparse(Aineq_i, Aineq_j, Aineq_v, nineqs, nvars);
P.CPLEX.Aeq = sparse(Aeq_i, Aeq_j, Aeq_v, neqs, nvars);
P.CPLEX.bineq=bineq;
P.CPLEX.beq=beq;
%objective
for j = 1:npossflows
    P.CPLEX.c(nass+nass*npossflows+j) = 1;
end

%variable bounds
%P.CPLEX.lb = -Inf*ones(numAss*numBatches+6,1);
%P.CPLEX.ub = Inf*ones(numAss*numBatches+6,1);

%print general problem parameters
fprintf('\t\tnumber of constraints = %i\n', neqs+nineqs);
fprintf('\t\t             equality = %i\n', neqs);
fprintf('\t\t           inequality = %i\n', nineqs);
fprintf('\t\tnumber of variables = %i\n', nvars);
fprintf('\t\t            integer = %i\n', sum(P.CPLEX.ctype == 'I'));
fprintf('\t\t             binary = %i\n', sum(P.CPLEX.ctype == 'B'));
fprintf('\t\t         continuous = %i\n', sum(P.CPLEX.ctype == 'C'));
