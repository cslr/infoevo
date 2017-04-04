
printf("%s\n", datestr(clock));

HISTORYSIZE = 30;
EPOCHS = 500; % generations

historyf1 = zeros(HISTORYSIZE,EPOCHS);
historyf2 = zeros(HISTORYSIZE,EPOCHS);

for h=1:HISTORYSIZE


% generates random genomes
M = 100; % population size
N = 20; % number of genes
G = rand(M,N) >= 0.45;

% number of generations
for g=1:EPOCHS

  printf("GENERATION %d/%d (%d/%d)\n", g, EPOCHS, h, HISTORYSIZE);
  fflush(stdout);
  
  t0 = time();
  
  Gnew = zeros(size(G));
  
  for k=1:M
    l1 = randi(M);
    l2 = randi(M);
    
    % 1. samples new genome
    gn = zeros(1, N);
    for n=1:N
      if(randi(2) == 1)
        gn(n) = G(l1, n);
      else
        gn(n) = G(l2, n);
      end
    end
    
    % fitness value
    
    p = sum(gn)/length(gn);
    
    % pa = sum(gn);
    % pb = length(gn) - pa;
    
    % f = binornd(1, p);
    % f = p;
    
    % acceptance function
    epsilon = 10e-10;
    
    q  = sum(G(k, :))/length(G(k, :)) + epsilon;
    
    % qa = sum(G(k,:));
    % qb = length(G(k,:)) - qa;
    % 
    % a = betapdf(f, pa, pb)/(betapdf(f, qa, qb) + epsilon);
    
    % a = binopdf(f, N, p)/(binopdf(f, N, q) + epsilon);
    
    % a = ((p/q)**(f)) * (((1-p)/(1-q))**(1 - f));
    
    if(p >= 0.5)
      p = p;
    else
      p = 1.0 - 0.95*p;
    end

    if(q >= 0.5)
      q = q;
    else
      q = 1.0 - 0.95*q;
    end


    
    a = (p/q); % **N (f=N)
    
    % printf("A = %f f = %f p = %f q = %f\n", a, f, p ,q);
    c = 1;
    a = a**c;
    
    if(a >= 1)
      Gnew(k, :) = gn;
      % printf("BETTER: %f\n", a);
      % fflush(stdout);
    else
      if(rand < a)
        Gnew(k,:) = gn;
        % printf("WORSE : %f\n", a);
        % fflush(stdout);
      else
        Gnew(k,:) = G(k,:);
        % printf("KEEP  : %f\n", a);
        % fflush(stdout);
      end   
    end
    
    
  end
  
  G = Gnew;
  
  % calculates average fitness value and best fitness value in population
  
  bestf = 0;
  averagef = 0;
  
  for k=1:M
    
    f = sum(G(k,:));
    
    if(f > bestf)
      bestf = f;
    end
    
    averagef = averagef + f/M;
    
  end
  
  historyf1(h, g) = bestf;
  historyf2(h, g) = averagef;
  
  t1 = time();
  
  minutesnow = (t1-t0)/60.0;
  hoursnow = (t1-t0)/3600.0;
  
  printf("Total time left: %f hours\n", hoursnow*EPOCHS*(HISTORYSIZE-h) + hoursnow*(EPOCHS-g));
  printf("Average fitness: %f\n", averagef);
  fflush(stdout);
  
  save -binary e3_history_data1
  
endfor

end

figure(1);
hold off;
plot(historyf1(1,:));
hold on;
for h=2:HISTORYSIZE
  plot(historyf1(h,:));
end

figure(2);
hold off;
plot(historyf2(1,:));
hold on;
for h=2:HISTORYSIZE
  plot(historyf2(h,:));
end

