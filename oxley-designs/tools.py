from fontTools.misc.cython import returns


class RefinementTask(object):
    def __init__(self, resolution, isinterface=False, newtagname=None):
        """
        interface: True if only the interface is to be refined
        tagname: the tag name to be used for elements below/inside interface
        resolution: target resolution
        """
        assert resolution >0
        self.tagname = tagname
        self.oppositetagid = None
        self.isinterface = isinterface
        self.resolution = resolution

    def isInterface(self):
        return self.isinterface
    def getTagName(self):
        return self.tagname
    def setTagId(self, tagid):
        self.tagid = tagid
    def getTagId(self):
        return self.tagid
    def setOppositeTageId(self, tagid):
        self.oppositetagid = tagid
    def getOppositeTageId(self):
        return self.oppositetagid
    def getResolution(self):
        return self.resolution
    def __call__(self, x):
        return self.check()

    def check(self, x):
        """
        return True if the x is inside/below the interface.
        """
        pass

class Sphere(RefinementTask):
    def __init__(self, center=(0.,0, 0.), radius=1., tagname="Sphere", resolution=1):
        super().__init__(resolution, isinterface=False, tagname=tagname)
        self.center = center
        self.radius = radius
    def check(self, x):
        d = math.sqrt((x[0]-self.center[0])**2 + (x[1]-self.center[1])**2+(x[2]-self.center[2])**2)
        return d <= self.radius
class PlaneInterface(RefinementTask):
    def __init__(self, origin=(0.,0, 400), normal=(0,0,1), tagname="Deep", resolution=10)
        super().__init__(resolution, isinterface=True, tagname=tagname)
        self.offset = inner(normal, origin)
        self.normal = normal
    def check(self, x):
        return inner(normal, x) <= -self.offset


class Refiner(object):
    def __init__(self, levels_max = 5):
        """
        levels_max: maximum number of refinemnet levels
        """
        assert levels_max > 0
        self.levels_max = levels_max
        self.tasks = []
    def add(self, tasks):
        """
        Add refinement tasks
        """
        if isinstance(tasks,list):
            self.tasks.extend(tasks)
        else:
            self.tasks.append(tasks)
    def getTasks(self):
        return self.tasks

    def refine(self, domain):
        tags = [ t.getTagName() for t in self.tasks if t.getTagName() is not None ]
        tagmap = domain.createTagmap(newtags=tags) -> tagmap[tagname]=tagid (chack for existing tags!)
        [ t.setTagID(tagmap(t.getTagName()) for t in self.tasks if t.getTagName() is not None]

        step = 0
        elements0 = domain.elements

        while step < self.levels_max:
            new_elements = []
            have_refined = False
            for e in elements0: iterate over tree (openmp?)
                splited_e = split(e)
                tag_e = e.getTagID()
                refine = False
                for t in self.tasks:
                    if e.size <= t.getResolution():
                        status = [ (e2, t.check(e2.x))  for e2 in splited_e ]
                        # subelements get tag according to side
                        if t.getTagName() is not None:
                            [e2.setTag(tagmap[t.getTagName()]) for e2, s in status if s ]
                            # if all subelements are inside, change tag of parent:
                            if all([s[1] for s in status]):
                                tag_e = tagmap[t.getTagName()]
                        # if we not deal with an interface then element is split if all
                        # subelements are on the blow the interface:
                        if  all([s[1] for s in status]):
                            if not t.isInterface():
                                refine = True
                        # if there is at least one subelement is below interface we split:
                        elif any([s[1] for s in status]):
                            refine = True
                if refine:
                    have_refined = True
                    new_elements.extend([s[1] for s in status])
                else:
                    new_elements.append(e.copy().setTagID(tag_e))
            -> do we need to iterate over surface elements??? tags are not updated!!!!

            if have_refined: ->MPI!!!
                - > update p4tree from new_elements
                elements0=new_elements
            else:
                break
            step+=1

    create domain and return







refner.add(Sphere(center=(0.,0, -400), radius=200, tagname="Anomaly", resolution=5))
refiner.add(PlaneInterface(origin=(0.,0, 400), normal=(0,0,1), tagname="Deep", resolution=10))
domain2=refiner(domain1)