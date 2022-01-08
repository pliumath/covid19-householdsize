function M = findHouseholdMembers(a,A)
%Find household members of an individual.

nummem = A(a,4);
indmem = A(a,3);
        
M = A((a-(indmem-1)):(a+(nummem-indmem)),:);

end

