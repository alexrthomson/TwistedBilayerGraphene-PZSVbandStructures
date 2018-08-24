function SixBand(V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectrum of 5-band honeycomb lattice model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin<1
		V = 0.0 ;
	end

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
	
	w = exp(2*pi*i/3) ; cw = conj(w) ;
	
	Vivc = 7*V ;
	
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
% pauli matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	sx = [ 0 1 ; 1 0 ];
	sy = [ 0 -i ; i 0 ];
	sz = [ 1 0 ; 0 -1 ];
	id = [ 1 0 ; 0 1 ];	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for tightbinding Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function [ham1,ham2]=symHam(p1,p2)
	
		phi11 = exp(i*(p1+p2)) ; cphi11 = conj(phi11) ;
		phi10 = exp(i*p1) ; cphi10 = conj(phi10) ;
		phi01 = exp(i*p2) ; cphi01 = conj(phi01) ; 
	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 1st valley
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% hopping between triangular lattice pz orbitals
		Hpz = tpz*( phi01 + phi11 + phi10 + cphi01 + cphi11 + cphi10  ) ;

		% hopping between triangular lattice p_+, p_- orbitals
		Cpm = tpmP*( phi01 + w*cphi11 + cw*phi10 ) + tpmM*( cphi01 + w*phi11 + cw*cphi10)  ;
		gpm = phi01 + phi11 + phi10 + cphi01 + cphi11 + cphi10 ;
		Hpm = tpm0*gpm.*id + [ 0, conj(Cpm) ; Cpm, 0 ] ;
	
		% hopping between kagome lattice s orbitals
		gk1 = [ 0, cphi10, 1.0 ; 1.0, 0.0, cphi01 ; phi11, 1.0, 0.0 ] ;
		gk2 = [ 0.0, cphi11, cphi10; cphi01, 0.0, phi10 ; phi01, phi11, 0.0 ] ;
		Hk = tk1.*( gk1 + gk1' ) + tk2.*( gk2 + gk2' ); 
	
		% hopping between p_z and p_{+/-} orbitals on triangular lattice 
		gpmzP = [ phi01 + w*cphi11 + cw*phi10 ; -( cphi01 + cw*phi11 + w*cphi10 ) ] ;
		gpmzM = [ cphi01 + w*phi11 + cw*cphi10 ; -( phi01 + cw*cphi11 + w*phi10 ) ] ;
		Cpmz = i*tpmzP.*gpmzP - i*tpmzM.*gpmzM ;
	
		% hopping between p_+/- orbtials on triangular lattice and the kagome s orbitals
		gkpmP = [ cphi10, cphi11 ; cw*cphi11, w ; w, cw*cphi10 ] ;
		gkpmM = [ cphi11, cphi10 ; cw, w*cphi11; w*cphi10, cw ] ;
		Ckpm = tkpmP.*gkpmP - tkpmM.*gkpmM ;
	
		% define Hamiltonian
		ham1 = zeros(6) ;
		ham1(1,1) = Hpz + muz ;
	
		ham1(2:3,2:3) = Hpm + mupm.*id ;
	
		ham1(4:6,4:6) = Hk + muk.*eye(3) ;
	
		ham1(2:3,1) = Cpmz ; 
		ham1(1,2:3) = Cpmz' ;
	
		ham1(4:6,2:3) = Ckpm ;
		ham1(2:3,4:6) = Ckpm' ;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% 2nd valley
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% hopping between triangular lattice pz orbitals
		HpzV2 = Hpz ;
	
		% hopping between triangular lattice p_+, p_- orbitals
		CpmV2 = tpmP*( phi01 + cw*cphi11 + w*phi10 ) + tpmM*( cphi01 + cw*phi11 + w*cphi10)  ;
		gpmV2 = gpm ;
		HpmV2 = tpm0*gpmV2.*id + [ 0, CpmV2 ; conj(CpmV2), 0 ] ;
		
		% hopping between kagome lattice s orbitals
		gk1V2 = [ 0, cphi10, 1.0 ; 1.0, 0.0, cphi01 ; phi11, 1.0, 0.0 ] ;
		gk2V2 = [ 0.0, cphi11, cphi10; cphi01, 0.0, phi10 ; phi01, phi11, 0.0 ] ;
		HkV2 = tk1.*( gk1V2 + gk1V2' ) + tk2.*( gk2V2 + gk2V2' ); 
		
		% hopping between p_z and p_{+/-} orbitals on triangular lattice 
		gpmzPV2 = [ -( cphi01 + w*phi11 + cw*cphi10 ) ; phi01 + cw*cphi11 + w*phi10 ] ;
		gpmzMV2 = [ -( phi01 + w*cphi11 + cw*phi10 ) ; cphi01 + cw*phi11 + w*cphi10 ] ;
		CpmzV2 = i*tpmzP.*gpmzPV2 - i*tpmzM.*gpmzMV2 ;
		
		% hopping between p_+/- orbtials on triangular lattice and the kagome s orbitals
		gkpmPV2 = [ cphi11, cphi10 ; cw, w*cphi11 ; w*cphi10, cw ] ;
		gkpmMV2 = [ cphi10, cphi11 ; cw*cphi11, w ; w, cw*cphi10 ] ;
		CkpmV2 = -tkpmP.*gkpmPV2 + tkpmM.*gkpmMV2 ;

		% define Hamiltonian
		ham2 = zeros(6) ;
		ham2(1,1) = HpzV2 + muz ;
	
		ham2(2:3,2:3) = HpmV2 + mupm.*id ;
	
		ham2(4:6,4:6) = HkV2 + muk.*eye(3) ;
	
		ham2(2:3,1) = CpmzV2 ; 
		ham2(1,2:3) = CpmzV2' ;
	
		ham2(4:6,2:3) = CkpmV2 ;
		ham2(2:3,4:6) = CkpmV2' ;
		

	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	EnKp = zeros(Nvec,Nvec,12) ;
	% EnKpV1 = zeros(Nvec,Nvec,6); EnKpV2 = EnKpV1 ;

	for n=1:Ntot
		
		kx = vkx(n); ky = vky(n) ;
		k1 = sqrt(3.0).*kx/2.0 - 0.5.*ky ;
		k2 = ky ;
		
		phi11 = exp(i*(k1+k2)) ;
		phi10 = exp(i*k1) ;
		phi01 = exp(i*k2) ;
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% tight binding hamiltonian 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
	
		[hamV1,hamV2] = symHam(k1,k2) ;
				
		ham0 = [ hamV1, zeros(6) ; zeros(6), hamV2 ] ;
	
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% IVC order
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
	
		hamIVC = Vivc.*kron( sx, eye(6) ) ;
		
		%*******************************************************
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% total hamiltonian
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%*******************************************************
	
		ham = ham0 + hamIVC ;
	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% diagonlize hamiltonian
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		[ U, En ] = eig(ham) ;
		dEn = diag(real(En));
		if ~issorted(dEn)
			[ dEn, I ] = sort(dEn) ;
			U = U(:,I) ;
		end
		% [ U, EnV1 ] = eig(hamV1) ;
		% dEnV1 = diag(real(EnV1));
		% [ dEnV1, I ] = sort(dEnV1,'descend') ;
		% U = U(:,I) ;
		
		% [ U, EnV2 ] = eig(hamV2) ;
		% dEnV2 = diag(real(EnV2));
		% [ dEnV2, I ] = sort(dEnV2,'descend') ;
		% % U = U(:,I) ;
		
		for j=1:6
			j1 = Nvec*Nvec*(j-1) + n ;
			EnKp(j1) = dEn(13-j) ;
		end
	
	end
	
%**************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find mins, maxes, and gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************************************************		

	for k=1:6
		Emins(k) = min(min(EnKp(:,:,k))) ;
		Emaxs(k) = max(max(EnKp(:,:,k))) ;
		
		% EminsV2(k) = min(min(EnKpV2(:,:,k))) ;
		% EmaxsV2(k) = max(max(EnKpV2(:,:,k))) ;
	end

	for k=1:5
		gaps(k) = Emins(k) - Emaxs(k+1) ;
		% gapsV2(k) = EminsV2(k) - EmaxsV2(k+1) ;
	end

	display(gaps) ;
	
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
	for k=1:12
		surf(vkx,vky,EnKp(:,:,k)) ;
	end
	xlabel( '$$k_x$$' );
	ylabel( '$$k_x$$' );
	view([1,1,1]);
	colorbar;
	title('All bands','Interpreter','latex');
	set( gcf, 'color', 'w' );
	
	

	f = figure;
	hold on ;
	for k=1:4
		surf(vkx,vky,EnKp(:,:,k)) ;
	end
	xlabel( '$$k_x$$' );
	ylabel( '$$k_x$$' );
	view([1,1,1]);
	colorbar;
	title('Flat bands close to CNP','Interpreter','latex');
	set( gcf, 'color', 'w' );


	% f = figure;
% 	hold on ;
% 	imagesc(pvec,pvec,EnKp(:,:,1))
% 	xlabel( '$$k_x$$' );
% 	ylabel( '$$k_y$$' );
% 	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
% 	colorbar;
% 	title('6-band: $$E_5$$','Interpreter','latex');
% 	set( gcf, 'color', 'w' );
% 	title('Top K valley bands','Interpreter','latex');
% 	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);
%
% 	f = figure;
% 	hold on ;
% 	imagesc(pvec,pvec,EnKpV2(:,:,1))
% 	xlabel( '$$k_x$$' );
% 	ylabel( '$$k_y$$' );
% 	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
% 	colorbar;
% 	title('6-band: $$E_5$$','Interpreter','latex');
% 	set( gcf, 'color', 'w' );
% 	title('Top K$$^\prime$$ valley bands','Interpreter','latex');
%
% 	f = figure;
% 	hold on ;
% 	imagesc(pvec,pvec,EnKpV1(:,:,2))
% 	xlabel( '$$k_x$$' );
% 	ylabel( '$$k_y$$' );
% 	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
% 	colorbar;
% 	title('6-band: $$E_5$$','Interpreter','latex');
% 	set( gcf, 'color', 'w' );
% 	title('Second K valley bands','Interpreter','latex');
% 	% set( findall(gca, 'Type', 'Line'), 'Linewidth', 1.5);
%
% 	f = figure;
% 	hold on ;
% 	imagesc(pvec,pvec,EnKpV2(:,:,2))
% 	xlabel( '$$k_x$$' );
% 	ylabel( '$$k_y$$' );
% 	line(bzX,bzY,'color','r', 'Linewidth', 2.0);
% 	colorbar;
% 	title('6-band: $$E_5$$','Interpreter','latex');
% 	set( gcf, 'color', 'w' );
% 	title('Second K$$^\prime$$ valley bands','Interpreter','latex');

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

