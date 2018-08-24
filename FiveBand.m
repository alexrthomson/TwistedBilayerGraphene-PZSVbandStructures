function [vkx,vky,E4,E5]=FiveBand
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	w = exp(2*pi*i/3) ;
	
	a = 0.25 ;
	b = 0.2 ;
	c = 0.1 ;
	d = 0.67 ;
	t0 = 80 ;

	muPz = -0.043.*t0 ;
	muPpm = 0.0 ;
	muEta = 0.05*t0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	E1 = zeros(Nvec) ; E2 = E1 ; E3 = E1 ; E4 = E1 ; E5 = E1 ;

	for n=1:Ntot
		
		kx = vkx(n); ky = vky(n) ;
		k1 = sqrt(3.0).*kx/2.0 - 0.5.*ky ;
		k2 = ky ;
		
		phi11 = exp(i*(k1+k2)) ;
		phi10 = exp(i*k1) ;
		phi01b = exp(-i*k2) ;
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% tight binding hamiltonian 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
	
		
		rho1 = [ i*a*( phi11 - phi10 ) ; b*phi11 + c*phi10 ; c*phi11 + b*phi10 ; ... 
			conj(d)*phi10 ; d ] ;
		rho2 = [ i*a*( 1.0 - phi11 ) ; w*( b + c*phi11 ) ; conj(w)*( c + b*phi11 ) ; ...
			conj(d) ; d ] ;
		rho3 = [ i*a*( phi10 - 1.0 ) ; conj(w)*( b*phi10 + c ) ; w*( c*phi10 + b ) ; ...
			conj(d)*phi01b ; d ] ;
			
		rho = [ rho1, rho2, rho3 ] ;
		
		ham = -t0.*rho*(rho') + diag( [ muPz, muPpm, muPpm, muEta, muEta ] ) ;
		
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
		E5(n) = dEn(5) ;

	
	end
	
	min1 = min(E1(:)) ; max1 = max(E1(:)) ;
	min2 = min(E2(:)) ; max2 = max(E2(:)) ;
	min3 = min(E3(:)) ; max3 = max(E3(:)) ;
	min4 = min(E4(:)) ; max4 = max(E4(:)) ;
	min5 = min(E5(:)) ; max5 = max(E5(:)) ;
	
	gaps = [ min2-max1, min3-max2, min4-max3, min5-max4 ];
	
	Emins = [min1,min2,min3,min4,min5]; 
	Emaxs = [max1,max2,max3,max4,max5]; 


	
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
	

	f = figure;
	hold on ;
	surf(vkx,vky,E1);
	surf(vkx,vky,E2);
	surf(vkx,vky,E3);
	surf(vkx,vky,E4);
	surf(vkx,vky,E5);
	xlabel( '$$k_x$$' );
	ylabel( '$$k_x$$' );
	view([1,1,1]);
	colorbar;
	% title('$$E_1$$','Interpreter','latex');
	set( gcf, 'color', 'w' );
	

	f = figure;
	hold on ;
	surf(vkx,vky,E4);
	surf(vkx,vky,E5);
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
	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
	colorbar;
	title('$$E_4$$','Interpreter','latex');
	set( gcf, 'color', 'w' );
	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);

	f = figure;
	hold on ;
	imagesc(pvec,pvec,E5)
	xlabel( '$$k_x$$' );
	ylabel( '$$k_y$$' );
	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
	colorbar;
	title('$$E_5$$','Interpreter','latex');
	set( gcf, 'color', 'w' );

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

