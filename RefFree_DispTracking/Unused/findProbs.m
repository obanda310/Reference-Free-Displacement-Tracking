function problems = findProbs(lub)

problems = find(lub(:,20)>0);
problems(:,2:6) = lub(problems,[15,17,19,22,20]);
