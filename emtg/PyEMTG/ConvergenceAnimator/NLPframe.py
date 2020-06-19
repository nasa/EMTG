#container and plotter for NLP movie frames
#holds vectors of Xrow and Frow objects

import Xrow
import Frow

class NLPframe(object):
    def __init__(self, frameFileName):
        self.parse_frame_file(frameFileName)

    def parse_frame_file(self, frameFileName):
        self.Xrows = []
        self.Frows = []

        passed_Fstart = False

        frameFileNameCell = frameFileName.split('_')
        self.frameIndex = int(frameFileNameCell[-1].replace('.csv',''))

        self.workdir = frameFileName.split('NLP')[0]

        frameFile = open(frameFileName, 'r')

        for line in frameFile:
            if 'Xindex' in line:
                continue
            elif 'Findex' in line:
                passed_Fstart = True
            else:
                if passed_Fstart:
                    self.Frows.append(Frow.Frow(line))
                else:
                    self.Xrows.append(Xrow.Xrow(line))

    def collect_violations(self):
        self.Xviolations = []
        self.XLogAbsViolations = []
        for Xrow in self.Xrows:
            self.Xviolations.append(Xrow.getViolation)
            self.XLogAbsViolations.append(Xrow.getLogAbsViolation)
                        
        self.Fviolations = []
        self.FLogAbsViolations = []
        for Frow in self.Frows:
            self.Fviolations.append(Frow.getViolation)
            self.FLogAbsViolations.append(Frow.getLogAbsViolation)

    def plot_violations(self):
        self.collect_violations()

        #and this is the part were Jacob needs to borrow pieces from PEATSA!!!

    def plot_F(self):
        #collect stuff
        self.Flowerbounds = []
        self.Fupperbounds = []
        self.Fvalues = []
        self.Fnormalized = []
        self.FLogAbsViolations = []
        self.Findex = []

        for Frow in self.Frows:
            self.Flowerbounds.append(Frow.getLowerBound())
            self.Fupperbounds.append(Frow.getUpperBound())
            self.Fvalues.append(Frow.getValue())
            if (self.Fupperbounds[-1] == self.Flowerbounds[-1]):
                self.Fnormalized.append(0.5)
            else:                
                self.Fnormalized.append((self.Fvalues[-1] - self.Flowerbounds[-1]) / (self.Fupperbounds[-1] - self.Flowerbounds[-1]))
            self.FLogAbsViolations.append(Frow.getLogAbsViolation())
            self.Findex.append(Frow.getIndex())

        #plot stuff
        
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        #make scatter plot
        plt.scatter(self.Findex[1:], self.Fnormalized[1:],c='blue')

        #save and clear scatter plot
        plt.savefig(self.workdir + 'F_frame' + str(self.frameIndex) + '.png',bbox_inches='tight')
        plt.clf()
        
        #make scatter plot
        plt.scatter(self.Findex[1:], self.FLogAbsViolations[1:],c='blue')

        #save and clear scatter plot
        plt.savefig(self.workdir + 'Flogabsviolation_frame' + str(self.frameIndex) + '.png',bbox_inches='tight')
        plt.clf()
        
    def plot_X(self):
        #collect stuff
        self.Xlowerbounds = []
        self.Xupperbounds = []
        self.Xvalues = []
        self.Xnormalized = []
        self.Xindex = []

        for Xrow in self.Xrows:
            self.Xlowerbounds.append(Xrow.getLowerBound())
            self.Xupperbounds.append(Xrow.getUpperBound())
            self.Xvalues.append(Xrow.getValue())
            if (self.Xupperbounds[-1] == self.Xlowerbounds[-1]):
                self.Xnormalized.append(0.5)
            else:
                self.Xnormalized.append((self.Xvalues[-1] - self.Xlowerbounds[-1]) / (self.Xupperbounds[-1] - self.Xlowerbounds[-1]))
            self.Xindex.append(Xrow.getIndex())

        #plot stuff
        
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        #make scatter plot
        plt.scatter(self.Xindex[1:], self.Xnormalized[1:],c='green')

        #save and clear scatter plot
        plt.savefig(self.workdir + 'X_frame' + str(self.frameIndex) + '.png',bbox_inches='tight')
        plt.clf()