function map = colormap_yyw(string)

%%%%string includes BlAqGrYeOrRe,BlGrYeOrReVi200,BlRe,BlueGreen14,BlueRed,BlWhRe
%%%%,example,GreenYellow,gui_default,rainbow,testcmap,ViBlGrWhYeOrRe

%%%%
pathdir = [string ,'.txt'];
color = load(pathdir);
color = color(3:end,:);
map = colormap(color);
