function [vkx,vky,E5,E6]=SixBandMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectrum of 5-band honeycomb lattice model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Nsite = 2000;
	
	% plim = 1.05*4.0*pi/3.0 ;
	% dp = 2*plim/(Nsite) ;
	% pvec=-plim:dp:plim ;
	% [vkx,vky]=meshgrid(pvec) ;
	%
	% Nvec = length(pvec) ;
	% Ntot = length(vkx(:)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	w = exp(2*pi*i/3) ;
	
	tk1 = 27 ;
	
	tpz = 0.17*tk1 ;
	tpm0 = -0.017*tk1 ;
	tpmP = -0.065*tk1 ;
	tpmM = -0.055*tk1 ;
	tk2 = 0.25*tk1 ;
	
	tpmzP = 0.095*tk1 ;
	tpmzM = 0.055*tk1 ;
	tkpmP = 0.6*tk1 ; 
	tkpmM = 0.2*tk1 ;
	
	muz = -6.0*tpz ;
	mupm = -0.23*tk1 + 3.*tpm0 ;
	muk = 0.25*tk1 - 4.0*(tk1+tk2) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make k-space line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%list points defining cut
	len = 4.0*pi/3.0 ;
	pts =	len.*[ 0, 1.0 ; 0.0, 0.0 ; sqrt(3.0)/4.0, 3.0/4.0 ; ] ;

	vecs = zeros(length(pts(:,1)),2) ;
	for k=1:length(pts(:,1))
	
		if k<length(pts(:,1)) 
			k2 = k+1 ;
		else
			k2 = 1 ;
		end
		
		legLens(k) = norm( pts(k2,:)-pts(k,:) )/len ;
		vecs(k,:) = ( pts(k2,:)-pts(k,:) )./norm( pts(k2,:)-pts(k,:) ) ;	
	end

	lenTot = sum(legLens) ;
	numPts = size(pts) ;
	numPts = numPts(1) ;

	kDist = 0.0 ;
	for k=1:numPts
		Nstep(k) = floor( legLens(k)*Nsite./lenTot ) ;
		dpx = vecs(k,1).*len*legLens(k)./Nstep(k) ;
		dpy = vecs(k,2).*len*legLens(k)./Nstep(k) ;

		vkxNow(1) = pts(k,1); vkyNow(1) = pts(k,2) ;
		for n=1:(Nstep(k)+1)
			vkxNow(n+1) = vkxNow(n) + dpx ;
			vkyNow(n+1) = vkyNow(n) + dpy ;
		end

		% vkxNow = (pts(k,1):dpx:pts(k+1,1)) ;
		%  	   vkyNow = (pts(k,2):dpy:pts(k+1,2)) ;

		if k==1
			vkx = vkxNow ; vky = vkyNow ;
		else
			vkx = [ vkx, vkxNow ] ;
			vky = [ vky, vkyNow ] ;
		end

		vkxNow = 0; vkyNow = 0 ;
	end

	Ntot = length(vkx(:)) ; 
	for k=1:Ntot

		if k==1 
			Kdist(k) = 0.0 ;
		else
			Kdist(k) = Kdist(k-1) + sqrt( (vkx(k)-vkx(k-1))^2.0 + (vky(k)-vky(k-1))^2.0 ) ;
		end

	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pauli matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	sx = [ 0 1 ; 1 0 ];
	sy = [ 0 -i ; i 0 ];
	sz = [ 1 0 ; 0 -1 ];
	id = [ 1 0 ; 0 1 ];	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	E1 = zeros(Ntot,1) ; E2 = E1 ; E3 = E1 ; E4 = E1 ; E5 = E1 ; E6 = E1 ;

	for n=1:Ntot
		
		kx = vkx(n); ky = vky(n) ;
		k1 = sqrt(3.0).*kx/2.0 - 0.5.*ky ;
		k2 = ky ;
		
		phi11 = exp(i*(k1+k2)) ;
		phi10 = exp(i*k1) ;
		phi01 = exp(i*k2) ;
		phi01b = exp(-i*k2) ;
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% tight binding hamiltonian 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
	
		gpz =  phi01 + phi11 + phi10  ; 
		Hpz = tpz*( gpz + conj(gpz) ) ;
	
		Cpm = tpmP*( phi01 + w*conj(phi11) + conj(w)*phi10 ) + tpmM*( conj(phi01) + w*phi11 + conj(w*phi10) ) ;
		gpm = phi01 + phi11 + phi10 ;
		Hpm = tpm0*( gpm + conj(gpm) ).*id + [ 0, conj(Cpm) ; Cpm, 0 ] ;
		
		gk1 = [ 0.0, conj(phi10), 1.0 ; 1.0, 0.0, conj(phi01) ; phi11, 1.0, 0.0 ] ;
		gk2 = [ 0.0, conj(phi11), conj(phi10); conj(phi01), 0.0, phi10 ; phi01, phi11, 0.0 ] ;
		Hk = tk1.*( gk1 + gk1' ) + tk2.*( gk2 + gk2' ); 
		
		gpmzP = [ phi01 + w*conj(phi11) + conj(w)*phi10 ; -( conj(phi01) + conj(w)*phi11 + w*conj(phi10) ) ] ;
		gpmzM = [ conj(phi01) + w*phi11 + conj(w*phi10) ; -( phi01 + conj(w*phi11) + w*phi10 ) ] ;
		Cpmz = i*tpmzP.*gpmzP - i*tpmzM.*gpmzM ;
		
		gkpmP = [ conj(phi10), conj(phi11) ; conj(w*phi11), w ; w, conj(w*phi10) ] ;
		gkpmM = [ conj(phi11), conj(phi10) ; conj(w), w*conj(phi11); w*conj(phi10), conj(w) ] ;
		Ckpm = tkpmP.*gkpmP - tkpmM.*gkpmM ;
		
		ham = zeros(6) ;
		ham(1,1) = Hpz + muz ;
		
		ham(2:3,2:3) = Hpm + mupm.*id ;
		
		ham(4:6,4:6) = Hk + muk.*eye(3) ;
		
		ham(2:3,1) = Cpmz ; 
		ham(1,2:3) = Cpmz' ;
		
		ham(4:6,2:3) = Ckpm ;
		ham(2:3,4:6) = Ckpm' ;
	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% diagonlize hamiltonian
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		[ U, En ] = eig(ham) ;
		dEn = diag(real(En));
		if ~issorted(dEn)
			[ dEn, I ] = sort(dEn) ;
			U = U(:,I) ;
		end
		E1(n) = dEn(1) ; E2(n) = dEn(2);
		E3(n) = dEn(3) ; E4(n) = dEn(4);
		E5(n) = dEn(5) ; E6(n) = dEn(6) ;

	
	end
	
	min1 = min(E1(:)) ; max1 = max(E1(:)) ;
	min2 = min(E2(:)) ; max2 = max(E2(:)) ;
	min3 = min(E3(:)) ; max3 = max(E3(:)) ;
	min4 = min(E4(:)) ; max4 = max(E4(:)) ;
	min5 = min(E5(:)) ; max5 = max(E5(:)) ;
	min6 = min(E6(:)) ; max6 = max(E6(:)) ;
	
	gaps = [ min2-max1, min3-max2, min4-max3, min5-max4, min6-max5 ];
	
	Emins = [min1,min2,min3,min4,min5,min6]; 
	Emaxs = [max1,max2,max3,max4,max5,max6]; 


	
%**************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************************************************	


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% fonts and font sizes
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	axFtSz = 18 ; labFtSz = 22 ;
	set(0,'defaulttextinterpreter','latex');
	set(0,'DefaultAxesFontName', 'CMU Serif');
	set(0,'defaultAxesFontSize',axFtSz);
	set(0,'defaultTextFontSize',labFtSz) ;
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% my colour scheme definitions
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	myOrange = [1,0.715,0];
	myGreen=[ 0.295, 0.8, 0.287 ];
	myRed = [ 1, 0.325, 0.407 ];
	myNavy = [ 0, 0.2, 0.4 ];
	myBlue = [0.6, 0.8, 1 ];
	

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% vertical line coords
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%vertical lines
	buff1 = (max6-min1)/20 ;
	buff2 = (max6-min5)/20 ;

	% points on large BZ
	vertX1 = len.*[ legLens(1), legLens(1) ];
	vertX2 = vertX1 + len.*[ legLens(2), legLens(2) ] ;
	vertX3 = vertX2 + len.*[ legLens(3), legLens(3) ] ;

	% length of line
	vertY1 = [min1-buff1,max6+buff1];
	vertY2 = [min5-buff2,max6+buff2];

	XticVec = [ 0, vertX1(1), vertX2(1), vertX3(1) ] ;

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	f = figure;
	hold on ;
	line(vertX1,vertY1,'LineStyle','--','Color',[0,0,0]+0.5 );
	line(vertX2,vertY1,'LineStyle','--','Color',[0,0,0]+0.5 );
	line(vertX3,vertY1,'LineStyle','--','Color',[0,0,0]+0.5 );
	plot(Kdist,E1,'Color',myNavy, 'Linewidth', 1.5 ) ;
	plot(Kdist,E2,'Color',myNavy, 'Linewidth', 1.5 ) ;
	plot(Kdist,E3,'Color',myNavy, 'Linewidth', 1.5 ) ;
	plot(Kdist,E4,'Color',myNavy, 'Linewidth', 1.5 ) ;	
	plot(Kdist,E5,'Color',myOrange, 'Linewidth', 1.5 ) ;	
	plot(Kdist,E6,'Color',myRed, 'Linewidth', 1.5 ) ;	
	xlabel( '$$k_x$$' );
	ylabel( '$$E$$ (meV)' );
	axis( [ 0, 2*pi*(3.0+sqrt(3.0))/3.0, min1-buff1, max6+buff1 ])
	xticks( XticVec ) ;
	% xticklabels( {'$$\Gamma$$', 'M', 'K', '$$\Gamma$$' });
	xticklabels( { 'K', '$$\Gamma$$', 'M', 'K' }) ;
	set( gca, 'TickLabelInterpreter', 'latex' );
	set( gcf, 'color', 'w' );
	
	
	f = figure;
	hold on ;
	line(vertX1,vertY2,'LineStyle','--','Color',[0,0,0]+0.5 );
	line(vertX2,vertY2,'LineStyle','--','Color',[0,0,0]+0.5 );
	line(vertX3,vertY2,'LineStyle','--','Color',[0,0,0]+0.5 );
	plot(Kdist,E5,'Color',myRed, 'Linewidth', 1.5 ) ;	
	plot(Kdist,E6,'Color',myBlue, 'Linewidth', 1.5 ) ;	
	xlabel( '$$k_x$$' );
	ylabel( '$$E$$ (meV)' );
	axis( [ 0, 2*pi*(3.0+sqrt(3.0))/3.0, min5-buff2, max6+buff2 ])
	xticks( XticVec ) ;
	% xticklabels( {'$$\Gamma$$', 'M', 'K', '$$\Gamma$$' });
	xticklabels( {'K', '$$\Gamma$$', 'M', 'K' }) ;
	set( gca, 'TickLabelInterpreter', 'latex' );
	set( gcf, 'color', 'w' );
	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);
	
	
end

