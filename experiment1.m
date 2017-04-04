
HISTORYSIZE = 30;
EPOCHS = 1000; % generations

historyf1 = zeros(HISTORYSIZE,EPOCHS);
historyf2 = zeros(HISTORYSIZE,EPOCHS);

for h=1:HISTORYSIZE


% generates random genomes
M = 1000; % population size
N = 20; % number of genes
G = rand(M,N) >= 0.5;

% number of generations that we go forward
for g=1:EPOCHS

  printf("GENERATION %d/%d (%d)\n", g, EPOCHS, h);
  fflush(stdout);
  
  t0 = time();
  
  for k=1:M
    l = randi(M);
    
    % 1. samples new genome
    gn = zeros(1, N);
    for n=1:N
      if(randi(2) == 1)
        gn(n) = G(k, n);
      else
        gn(n) = G(l, n);
      end
    end
    
    % fitness value
    
    p = sum(gn)/length(gn);
    
    %% f = binornd(N, p);
    
    f = binornd(N, p);
    
    % f = 1;
    
    % acceptance function
    epsilon = 10e-10;
    
    q  = sum(G(k, :))/length(G(k, :)) + epsilon;
    
    % a = binopdf(f, N, p)/(binopdf(f, N, q) + epsilon);
    
    a = ((p/q)**(f)) * (((1-p)/(1-q))**(N - f));
    
    %% a = (p/q)**f * ((1-p)/(1-q))**(1-f);
    
    c = 10;
    a = a**c;
    
    if(a >= 1)
      G(k, :) = gn;
      %printf("Selecting better solution: %f\n", a);
      %fflush(stdout);
    else
      if(rand < a)
        G(k,:) = gn;
        %printf("Selecting worse solution: %f\n", a);
        %fflush(stdout);
      end
       
    end
    
    
  end
  
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
  
  printf("Computation time: %f minutes\n", minutesnow);
  printf("Average fitness: %f\n", averagef);
  printf("Total time estimate: %f hours\n", hoursnow*EPOCHS*HISTORYSIZE);
  fflush(stdout);
  
  save -binary history_data10
  
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

