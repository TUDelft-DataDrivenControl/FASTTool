function writebladed(filename,u,v,w,x,y,z,U)
% syntax function writebladed(filename,u,v,w,x,y,z,U)

[Nx,Ny,Nz]=size(u);

WW='''w''';
eval(['fid=fopen(''',filename,'.wnd''',',',WW,');']);

New=-99;
fwrite(fid,New,'int16');

Spec=3;
fwrite(fid,Spec,'uint16');

delta_x=x(2)-x(1);
% BLADED: first element at right bottom corner, i.e. maximum y
delta_y=-(y(2)-y(1));
delta_z=z(2)-z(1);

fwrite(fid,delta_z,'float32');
fwrite(fid,delta_y,'float32');
fwrite(fid,delta_x,'float32');

Nx2=Nx/2;
Nx2=fwrite(fid,Nx2,'uint32');

fwrite(fid,U,'float32');

Lux=70.7;
Luy=Lux/2;
Luz=Lux/2;
fwrite(fid,Luz,'float32');
fwrite(fid,Luy,'float32');
fwrite(fid,Lux,'float32');

fwrite(fid,1,'uint32');
Seed=1;
fwrite(fid,Seed,'uint32');

fwrite(fid,Nz,'uint32');
fwrite(fid,Ny,'uint32');

% length scales
fwrite(fid,100,'float32');
fwrite(fid,100,'float32');
fwrite(fid,100,'float32');
fwrite(fid,100,'float32');
fwrite(fid,100,'float32');
fwrite(fid,100,'float32');

Nplane=Ny*Nz*3;
Ind1=1:3:Nplane;
Ind2=2:3:Nplane;
Ind3=3:3:Nplane;
% scaling
u=u*1000;
v=v*1000;
w=w*1000;
uvw=zeros(1,Nplane);
for i=1:Nx,
   uu=reshape(u(i,:,:),1,Ny*Nz);
   vv=reshape(v(i,:,:),1,Ny*Nz);
   ww=reshape(w(i,:,:),1,Ny*Nz);
   uvw(Ind1)=uu;
   uvw(Ind2)=vv;
   uvw(Ind3)=ww;
   fwrite(fid,uvw,'int16'); 
end

fclose(fid);
