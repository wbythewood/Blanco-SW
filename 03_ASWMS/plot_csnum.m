clear;

setup_parameters;
comp = parameters.component;
files = dir(['./03_ASWMS/eikonal/*',comp,'.mat']);


for ie = 1:length(files)
	disp(files(ie).name)
	filename = ['./03_ASWMS/eikonal/',files(ie).name];
	load(filename);
	good_bad_ratio(:,ie) = [eventphv.goodnum] ./ [eventphv.badnum];
end
good_bad_ratio(find(good_bad_ratio>2)) = 2;

figure(22)
clf
N=3; M = floor(size(good_bad_ratio,1)/N)+1;
for ip = 1:size(good_bad_ratio,1)
	subplot(N,M,ip)
	hist(good_bad_ratio(ip,:));
end

