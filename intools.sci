function [x,y] = icsvRead (name, separator, skipheader)
    xy = csvRead(name, separator, [], [], [], [], [], skipheader);
    x = xy(:,1);
    y = interval(xy(:,2),xy(:,3));
endfunction

function iplot (x,y)
    errbar(x,mid(y),rad(y),rad(y));
    plot(x,mid(y),'ko');
    plot(x,sup(y),'k--');
    plot(x,inf(y),'k--');
endfunction


// Builds interval regression LP problem
// Reference: 
//    Yu et al. / Fuzzy Sets nd Systems 105 (1999) 429-436,
//         Eq. (14)-(15), (16)-(19)
//    NB! There is an error in the paper: 
//        the sign of a_{0w*} in Eq. (17) must be minus!
function [c, A, b, ci, cs] = ibuildLP (x,y)
    [n q] = size(x);
    
    // Objective function's coefficients for variables 
    // [a0c a1c ... aqc a0w a1w ... aqw]
    //   aic = mid(A_i), i=0,...,q
    //   aiw = rad(A_i), i=0,...,q    
    //
    // [a0c a1c ... aqc a0w     a1w    ...    aqw    ]
    // [ 0   0  ...  0   n  sum(|x1j|) ... sum(|x1j|)]
    c = [zeros(1,q+1) n sum(abs(x),'r')];

    A = [-ones(n,1) -x  ones(n,1) abs(x);
          ones(n,1)  x  ones(n,1) abs(x)];
          
    b = [-inf(y); sup(y)];
    
    //Infty=number_properties('huge');
    Infty=100000;
    ci = [-Infty*ones(1,q+1) zeros(1,q+1)]';
    cs = Infty*ones(1,2*(q+1))';
endfunction


// Solves interval regression LP problem
// Reference: 
//    Yu et al. / Fuzzy Sets nd Systems 105 (1999) 429-436,
//         Eq. (14)-(15), (16)-(19)
function [a,exitcode] = isolveLP (c, A, b, ci, cs)
    //[x [,iact [,iter [,f]]]] = qpsolve(Q,p,C,b,ci,cs,me)
//    Qn = size(c,'c');
//    Q = zeros(Qn,Qn);
//    disp(Q)
//    [acw] = qpsolve(Q,c',A,b,ci,cs,0);
    
    //[xopt,fopt,exitflag] = karmarkar(Aeq,beq,c,x0,rtolf,gam,maxiter,outfun,A,b,lb,ub)
    //[acw,fopt,exitflag] = karmarkar([],[],-c',[],[],[],100,[],A,b,ci,cs)
    //disp(exitflag)
    
    a = zeros(c)
    exitcode = 1
    try
      [acw,lagr,flow]=linpro(-c',A,b,ci,cs,0)
    catch
      [error_message,error_number]=lasterror(%t)
      select error_number
          case 123 then
                exitcode = -1  
          case 127 then
                exitcode = -2
          else      
                exitcode = -3
                disp(error_number)
                disp(error_message)
      end 
      return 
    end

    //disp(acw)
    a = interval(acw(1:$/2)-acw($/2+1,$),acw(1:$/2)+acw($/2+1,$));
    //a=acw
    
endfunction
