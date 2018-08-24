function [gaps,Emins,Emaxs]=FiveBandMagOrd(M,Ma,chemPot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectrum of 5-band honeycomb lattice model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Nsite = 100;
	
	plim = 1.05*4.0*pi/3.0 ;
	dp = 2*plim/(Nsite) ;
	pvec=-plim:dp:plim ;
	[vkx,vky]=meshgrid(pvec) ;

	Nvec = length(pvec) ;
	Ntot = length(vkx(:)) ;
	
	Qyh = 2.0*pi/3.0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	w = exp(2*pi*i/3) ; cw = conj(w) ;
	
	t0 = 80 ;
	
	a = 0.25 ;
	b = 0.2 ;
	c = 0.1 ;
	d = 0.67 ; ccd = conj(d) ;


	muPz = -0.043.*t0 ;
	muPpm = 0.0 ;
	muEta = 0.05*t0 ;
	
	
	M0 = M*7 ; % 7 meV is the bandwidth these parameters give
	M1 = Ma*7 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pauli matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	sx = [ 0 1 ; 1 0 ];
	sy = [ 0 -i ; i 0 ];
	sz = [ 1 0 ; 0 -1 ];
	id = [ 1 0 ; 0 1 ];	
	
	sPio3 = cos(pi/3.0).*sx + sin(pi/3.0).*sy ;
	s5Pio6 = cos(5.0*pi/6.0).*sx + sin(5.0*pi/6.0).*sy ;
	
	% project onto triangular lattice sites
	Pt = diag([ 1 1 1 0 0]);
	
	Ph0 = diag([0,0,0,1,1]) ;
	Phz = diag([0,0,0,1,-1]) ;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for tightbinding Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
asdfasg

gfds
	% matrix to implement time reversal transformation
	mtrs = eye(5); 
	mtrs(2:3,2:3) = -sx ;
	
	function hamOut=symHam(p1,p2)
		
		phi11 = exp(i*(p1+p2)) ; cphi11 = conj(phi11) ;
		phi10 = exp(i*p1) ; cphi10 = conj(phi10) ;
		phi01 = exp(i*p2) ; cphi01 = conj(phi01) ;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 1st valley
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		rho1v1 = [ i*a*( phi11 - phi10 ) ; b*phi11 + c*phi10 ; c*phi11 + b*phi10 ; ... 
			ccd*phi10 ; d ] ;
		rho2v1 = [ i*a*( 1.0 - phi11 ) ; w*( b + c*phi11 ) ; conj(w)*( c + b*phi11 ) ; ...
			ccd ; d ] ;
		rho3v1 = [ i*a*( phi10 - 1.0 ) ; conj(w)*( b*phi10 + c ) ; w*( c*phi10 + b ) ; ...
			ccd*phi01b ; d ] ;
			
		rhov1 = [ rho1v1, rho2v1, rho3v1 ] ;
		
		ham1 = rhov1*(rhov1') ;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 2st valley
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		rho1v2 = [ -i*a*( phi11 - phi10 ) ; b*phi11 + c*phi10 ; c*phi11 + b*phi10 ; ... 
			d*phi10 ; ccd ] ;
		rho2v2 = [ -i*a*( 1.0 - phi11 ) ; cw*( b + c*phi11 ) ; w*( c + b*phi11 ) ; ...
			d ; ccd ] ;
		rho3v2 = [ -i*a*( phi10 - 1.0 ) ; c*( b*phi10 + c ) ; cw*( c*phi10 + b ) ; ...
			d*phi01b ; ccd ] ;
			
		rhov2 = [ rho1v2, rho2v2, rho3v2 ] ;
		
		ham2 = rhov2*(rhov2') ;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% together
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		hamOut = [ ham1, eye(5) ; eye(5), ham2 ] ;
	
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	EnKp = zeros(Nvec,Nvec,24) ;

	for n=1:Ntot
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% momenta
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		kx = vkx(n); ky = vky(n) ;

		k1 = sqrt(3.0).*kx/2.0 - 0.5.*ky ;
		k2 = ky ;

		
		k1P = sqrt(3.0).*kx/2.0 - 0.5.*(ky+Qyh) ;
		k2P = ky + Qyh ;
		k1M = sqrt(3.0).*kx/2.0 - 0.5.*(ky-Qyh) ;
		k2M = ky - Qyh ;
		
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% tight binding hamiltonian 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
	
		h0P = symHam(k1P,k2P) ; 
		h0M = symHam(k1M,k2M) ;
		
		ham0 = -t0.*[ h0P, zeros(5) ; zeros(5), h0M ] + kron(id, diag( [ muPz, muPpm, muPpm, muEta, muEta ] ) ) ;
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% magnetic hamiltonian
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
		
		hamMt = M0.*kron( sx, Pt ) ;
		
		hamH2m = 2.0.*M1.*( (cos(k1+k2)-0.5.*(cos(k1)-cos(k2))).*kron(sPio3,Ph0) + ...
		 	 (sqrt(3.0)/2.0).*(cos(k1)+cos(k2)).*kron(s5Pio6,Phz) ) ;
			 
	    hamM = hamMt + hamH2m ;
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% total hamiltonian
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
		
		ham = ham0 + hamM ;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% diagonlize hamiltonian
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		[ U, En ] = eig(ham) ;
		dEn = diag(real(En));
		[ dEn, I ] = sort(dEn,'descend') ;
		U = U(:,I) ;
		
		E1(n) = dEn(1) ; E2(n) = dEn(2);
		E3(n) = dEn(3) ; E4(n) = dEn(4);
		E5(n) = dEn(5) ;
		E6(n) = dEn(6) ; E7(n) = dEn(7);
		E8(n) = dEn(8) ; E9(n) = dEn(9);
		E10(n) = dEn(10) ;

	
	end
	
	min1 = min(E1(:)) ; max1 = max(E1(:)) ;
	min2 = min(E2(:)) ; max2 = max(E2(:)) ;
	min3 = min(E3(:)) ; max3 = max(E3(:)) ;
	min4 = min(E4(:)) ; max4 = max(E4(:)) ;
	min5 = min(E5(:)) ; max5 = max(E5(:)) ;
	min6 = min(E6(:)) ; max6 = max(E6(:)) ;
	min7 = min(E7(:)) ; max7 = max(E7(:)) ;
	min8 = min(E8(:)) ; max8 = max(E8(:)) ;
	min9 = min(E9(:)) ; max9 = max(E9(:)) ;
	min10 = min(E10(:)) ; max10 = max(E10(:)) ;
	
	gaps = [ min1-max2, min2-max3,min3-max4,min4-max5,min5-max6,min6-max7,min7-max8,min8-max9,min9-max10 ];

	display(gaps);
	
	Emins = [min1,min2,min3,min4,min5,min6,min7,min8,min9,min10]; 
	Emaxs = [max1,max2,max3,max4,max5,max6,max7,max8,max9,max10]; 


	
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
	% BZ plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bzX = (4*pi/3).*[ -sqrt(3)/2, 0. ; 0., sqrt(3)/2 ; sqrt(3)/2, sqrt(3)/2 ; ...
	sqrt(3)/2, 0. ; 0., -sqrt(3)/2; -sqrt(3)/2, -sqrt(3)/2 ];
	bzY = (4*pi/3).*[ 1/2, 1 ; 1., 1/2 ; 1/2, -1/2; ...
	-1/2, -1 ; -1, -1/2; -1/2, 1/2 ];

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% plot
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%
	f = figure;
	hold on ;
	surf(vkx,vky,E1);
	surf(vkx,vky,E2);
	surf(vkx,vky,E3);
	surf(vkx,vky,E4);
	surf(vkx,vky,E5);
	surf(vkx,vky,E6);
	surf(vkx,vky,E7);
	surf(vkx,vky,E8);
	surf(vkx,vky,E9);
	surf(vkx,vky,E10);
	xlabel( '$$k_x$$' );
	ylabel( '$$k_x$$' );
	view([1,1,1]);
	colorbar;
	% title('$$E_1$$','Interpreter','latex');
	set( gcf, 'color', 'w' );
	

	f = figure;
	hold on ;
	surf(vkx,vky,E1);
	surf(vkx,vky,E2);
	surf(vkx,vky,E3);
	surf(vkx,vky,E4);
	xlabel( '$$k_x$$' );
	ylabel( '$$k_x$$' );
	view([1,1,1]);
	colorbar;
	% title('$$E_1$$','Interpreter','latex');
	set( gcf, 'color', 'w' );

	f = figure;
	hold on ;
	imagesc(pvec,pvec,E4)
	xlabel( '$$k_x$$' );
	ylabel( '$$k_y$$' );
	line(bzX,bzY,'color','k', 'Linewidth', 2.0);
	hold on ;
	contour(vkx,vky,E4,[chemPot,chemPot], 'Linewidth', 2.0,'color','r') ;
	colorbar;
	title('$$E_4$$','Interpreter','latex');
	set( gcf, 'color', 'w' );
	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);

	f = figure;
	hold on ;
	imagesc(pvec,pvec,E3)
	xlabel( '$$k_x$$' );
	ylabel( '$$k_y$$' );
	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
	hold on ;
	contour(vkx,vky,E3,[chemPot,chemPot], 'Linewidth', 2.0,'color','r') ;
	colorbar;
	title('$$E_3$$','Interpreter','latex');
	set( gcf, 'color', 'w' );
	
	% f = figure;
% 	hold on ;
% 	imagesc(pvec,pvec,E9)
% 	xlabel( '$$k_x$$' );
% 	ylabel( '$$k_y$$' );
% 	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
% 	colorbar;
% 	title('$$E_9$$','Interpreter','latex');
% 	set( gcf, 'color', 'w' );
% 	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);
%
% 	f = figure;
% 	hold on ;
% 	imagesc(pvec,pvec,E10)
% 	xlabel( '$$k_x$$' );
% 	ylabel( '$$k_y$$' );
% 	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
% 	colorbar;
% 	title('$$E_{10}$$','Interpreter','latex');
% 	set( gcf, 'color', 'w' );

	%
	% f = figure;
	% hold on ;
	% imagesc(pvec,pvec,E3)
	% xlabel( '$$k_x$$' );
	% ylabel( '$$k_y$$' );
	% line(bzX,bzY,'color','r');
	% colorbar;
	% title('$$E_3$$','Interpreter','latex');
	% set( gcf, 'color', 'w' );
	%
	% f = figure;
	% hold on ;
	% imagesc(pvec,pvec,E4)
	% xlabel( '$$k_x$$' );
	% ylabel( '$$k_y$$' );
	% line(bzX,bzY,'color','r');
	% colorbar;
	% title('$$E_4$$','Interpreter','latex');
	% set( gcf, 'color', 'w' );
	
	
	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);
	%
	% f = figure;
	% hold on ;
	% imagesc(vk1,vk2,E4)
	% xlabel( '$$k_1$$' );
	% ylabel( '$$k_2$$' );
	% colorbar;
	% title(leg,'$$E_4$$','Interpreter','latex');
	% set( gcf, 'color', 'w' );
	% % set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);

end

