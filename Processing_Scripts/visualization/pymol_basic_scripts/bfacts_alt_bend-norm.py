from pymol import cmd, stored, math
	
def loadBfacts(structurename, source=pathbf+'/'+structurename+"_bend-norm.txt"):
    newb = [i for i in open(source).read().split()]
    for x in range(0, len(newb)):
        cmd.alter("%s and resi %s"%(structurename,x+1), "b=%s"%newb[x+1])
    return newb


def visualtest(structurename, bfactlist):
    cmd.show_as("cartoon",structurename)
    cmd.cartoon("putty", structurename)
    cmd.set("cartoon_putty_scale_min", min(bfactlist))
    cmd.set("cartoon_putty_scale_max", max(bfactlist))
    cmd.set("cartoon_putty_transform", 0)
    cmd.set("cartoon_putty_radius", 0.2)
    cmd.spectrum("b","rainbow", "%s" %structurename)
    cmd.recolor()
    return
    
