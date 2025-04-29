#!/usr/bin/env python3

import imas
import argparse

import os
import re

import logging
log = logging.getLogger(__name__)

class EQDSKIO(object):
    """EQDSK interface implementation

    The following quantities are read from a ``EQDSK G``.

    Attributes:
        PSIRZ (arr): Array of plasma psi values in units Web/rad.
        CURRENT (float): Plasma current in Ampere
        RDIM (float): Horizontal dimension in meter of computational box
        ZDIM (float): Vertical dimension in meter of computational box
        NW (int): Number of horizontal R grid points
        NH (int): Number of vertical Z grid points
        LIMITR (int): Number of limiter points
        RLIM (arr): R of surrounding limiter contour in meter
        ZLIM (arr): Z of surrounding limiter contour in meter
        NBBBS (int): Number of boundary points
        RBBBS (arr): R of boundary points in meter
        ZBBBS (arr): Z of boundary points in meter
        RMAXIS (arr): R of magnetic axis in meter
        ZMAXIS (arr): Z of magnetic axis in meter
        FPOL (arr): Poloidal current function in m-T, F=RB_T on flux grid
        PRES (arr): Plasma pressure in nt/m^2 on uniform flux grid
        FFPRIM (arr): FF'(Psi) in (mT)^2/(Weber/rad) on uniform flux grid
        PPRIME (arr): P'(Psi) in (nt/m^2) / (Weber/rad) on uniform flux grid
        QPSI (arr): q values on uniform flux grid from axis to boundary
        SIMAG (arr): poloidal flux at magnetic axis in Weber /rad
        SIBRY (arr): poloidal flux at the plasma boundary in Weber/rad
        RCENTR (float): R in meter of vacuum toroidal magnetic field BCENTR
        BCENTR (float): Vacuum toroidal magnetic field in Tesla at RCENTR
        RLEFT (float): Minimum R in meter of rectangular computational box
        ZMID (float): Z of center of computational box in meter

    They are all accessible with the following example

    ``` python
       x = EQDSKIO()
       eqdskFile = '/path/to/eqdsk/format/file'
       # If SMITER is located on your machine
       x.openAndRead(eqdskFile)

       # If SMITER is used remotely you have to manually pass the contents
       # of the EQDSK file. In this case:
       with open(eqdskFile) as f:
           data = f.read()
       ok = x.read(data)
       # Then set the name of the eqdsk file
       if ok:
           # Reading Successful
           eqdskName = eqdskFile.rpslit('/', 1)[-1] # Get the name of the file
           x.setName(eqdskFile)
       else:
           print('Trouble reading EQDSK file %s' % eqdskFile)

       # If the file was read without a problem, you can access, i.e, the
       # geometry of the limiter
       x.RLIM
       x.ZLIM
    ```
    """

    # Pattern for reading integer values separated with spaces
    integerPattern = r'(?:^|(?<=\s))-?[0-9]+(?:$|(?=\s))'
    intPatt = re.compile('.*' + integerPattern + '.*')
    # Pattern for reading scientific numbers
    scientificPattern = r'[+-]?[0-9]\.[0-9]{6,9}[eE][+-][0-9]{2}'

    pattern = u'(?:%s|%s)' % (integerPattern, scientificPattern)
    genPatt = re.compile(pattern)

    def __init__(self, file=''):
        self.name = ''
        self.logMessage = ''
        self.resetValues()

        if file:
            ok = self.openAndRead(file)
            if not ok:
                log.error(f'Failed to load file: {file}')

    def resetValues(self):
        self.INPUT_EQDSK = ''
        self.HEADER = ''

        self.PSIRZ = [] # Plasma psi WEBER

        self.CURRENT = None

        self.RDIM = None # Size dimension for R and Z
        self.ZDIM = None

        self.NW = None  # Dimension of grid
        self.NH = None

        self.LIMITR = None  # Number of limiter points
        self.RLIM = []  # rlim
        self.ZLIM = []  # zlim

        self.NBBBS = None  # Number of Boundary points
        self.RBBBS = []  # rbbbs
        self.ZBBBS = []  # zbbs

        self.RMAXIS = []  # Position of toroidal magnetic
        self.ZMAXIS = []  # field
        self.FPOL = []
        self.PRES = []
        self.FFPRIM = []
        self.PPRIME = []
        self.QPSI = []

        self.SIMAG = None  # Flux value at magnetix axis
        self.SIBRY = None  # Flux value on boundary

        self.RCENTR = None  # Reference value: position R, mag field B
        self.BCENTR = None  # B value in RCENTR

        self.RLEFT = None  # Minimum of rectangular computational box.
        self.ZMID = None  # Vertical dimension of computational box
        """
        Or in another case:

        The grid goes from :

        RLEFT to RLEFT + RDIM
        -ZDIM / 2  to  ZDIM / 2
        """

        self.eqdskString = ''
        self.successfullRead = False

    def getName(self) -> str:
        """Returns the name.

        Returns:
            name (str): String which has whites-paces replaced with `_`.
        """
        return self.name.replace(' ', '_')

    def setName(self, name: str) -> None:
        """Set a name for eqdsk.. It should be the same as the name displayed
        to the **Salome Object** in **Object browser**.

        Arguments:
            name (str): Name for the case.
        """
        self.name = name

    def openAndRead(self, filePath: str) -> bool:
        if not os.access(filePath, os.F_OK | os.R_OK):
            log.error(f'File {filePath} does not exist or lacking permission to read it!')
            return False

        data = ''
        with open(filePath, 'r') as f:
            data = f.read()
        # Setting name
        self.setName(os.path.basename(filePath))
        return self.read(data)

    def logValue(self, msg, value, t='f') -> None:
        if t == 'd':
            formatedMsg ='{0:>20} %d'.format(msg + ':') % value
        else:
            formatedMsg ='{0:>20} %f'.format(msg + ':') % value
        self.logMessage += formatedMsg + '\n'

    def log(self, msg: str) -> None:
        """Function used to save the parsing information of the EQDSK file.
        """
        self.logMessage += msg + '\n'

    def getLog(self) -> str:
        return self.logMessage

    def read(self, eqdskString: str) -> bool:
        """Passing string instead of a filename, to avoid any troubles if this
        is used on a server.

        This function parses EQDSK G format.
        """
        self.resetValues()
        self.eqdskString = eqdskString
        # Reset the log
        self.logMessage = 'Reading eqdsk file %s\n' % self.name

        LINES = eqdskString.splitlines()

        HEADER = LINES[0]
        self.HEADER = HEADER
        # regexHeader = '^.+[0-9] +?([0-9]+) +?([0-9]+)'
        regexHeader = r'\s[0-9]\s+([0-9]+)\s+([0-9]+)'

        result = re.findall(regexHeader, HEADER)
        # self.log("Reading 1st line of EQDSK.")
        if result:
            nw, nh = result[0]
            self.NW = int(nw)
            self.NH = int(nh)
            self.logValue('Width', self.NW, t='d')
            self.logValue('Height', self.NH, t='d')
        else:
            self.log("Couldn't read the dimensions from the 1st line.")
            self.log("LINE:\n%s" % HEADER)
            return False
        # self.log('Done reading 1st line.')
        # self.log("Reading 2nd line of EQDSK")
        i = 1
        result = self.readOneLine(LINES[i])

        if result == -1:
            # self.log("[INFO] 2nd line must be a comment")
            while 1:
                # self.log("[INFO] Filtering the comments.")
                i += 1

                result = self.readOneLine(LINES[i])
                if result != -1:
                    break

                # Precaution, if there are 100 wrong lines, just end the
                # reading
                if i > 100:
                    return False
        try:
            self.RDIM = result[0]
            self.ZDIM = result[1]
            self.RCENTR = result[2]
            self.RLEFT = result[3]
            self.ZMID = result[4]

            self.logValue('Dimension R', self.RDIM)
            self.logValue('Dimension Z', self.ZDIM)
            self.logValue('Center R', self.RCENTR)
            self.logValue('Left R', self.RLEFT)
            self.logValue('Middle Z', self.ZMID)
            # self.log('Done reading 2nd line.')
            # self.log("Reading 3rd line of EQDSK")
            i += 1
            result = self.readOneLine(LINES[i])
            self.RMAXIS = result[0]
            self.ZMAXIS = result[1]
            self.SIMAG = result[2]
            self.SIBRY = result[3]
            self.BCENTR = result[4]

            self.logValue('Magnetic axis position R',
                           self.RMAXIS)
            self.logValue('Magnetic axis position Z',
                           self.ZMAXIS)
            self.logValue('Flux at magnetic axis', self.SIMAG)
            self.logValue('Flux at boundary', self.SIBRY)
            self.logValue('Center B', self.BCENTR)
            # self.log("Done reading 3rd line.")

            i += 1
            # self.log('Reading 4th line of EQDSK')
            result = self.readOneLine(LINES[i])
            self.CURRENT = result[0]
            self.logValue('Current', self.CURRENT)

            msg = "The 2nd value on line 4 of the EQDSK is not equal to " \
                  "the flux at magnetic axis read from the 3rd line"
            assert self.SIMAG == result[1], msg
            DUMMY = result[2]
            msg = "The 4th value on line 4 of the EQDSK is not equal to the"\
                  " magnetic axis R position, read from the 3rd line!"
            assert self.RMAXIS == result[3], msg
            DUMMY = result[4]
            # self.log("Done reading 4th line.")
            i += 1
            # self.log("Reading 5th line of EQDSK")
            result = self.readOneLine(LINES[i])
            msg = "The 1st value on line 5 is not equal to the magnetic " \
                  " axis Z position, read from the 3rd line!"
            assert self.ZMAXIS == result[0], msg
            DUMMY = result[1]
            msg = " The 3rd value on line 5 is not equal to the flux at " \
                  " boundary, read from the 3rd line"
            assert self.SIBRY == result[2], msg
            DUMMY = result[3]
            DUMMY = result[4]
            # self.log("Done reading 5th line.")

            line_index = i + 1
            self.log("Reading values from 5th line on...")
            P1, P2 = self.readValues(LINES[line_index:])
            N_VALUES = len(P1) + len(P2)
            TOTAL_VALUES = 5 * self.NW + self.NH * self.NW + 2

            _I = 0
            self.FPOL = P1[_I: _I + self.NW]
            _I += self.NW

            self.PRES = P1[_I: _I + self.NW]
            _I += self.NW

            self.FFPRIM = P1[_I: _I + self.NW]
            _I += self.NW

            self.PPRIME = P1[_I: _I + self.NW]
            _I += self.NW

            self.PSIRZ = []
            for i in range(self.NH):
                self.PSIRZ.append(P1[_I: _I + self.NW])
                _I += self.NW

            self.QPSI = P1[_I: _I + self.NW]
            _I += self.NW

            _I = 0
            self.NBBBS = int(P2[_I])
            _I += 1

            self.LIMITR = int(P2[_I])
            _I += 1
            self.logValue('N boundary', self.NBBBS, t='d')
            self.logValue('N limiter', self.LIMITR, t='d')
            self.log("Read %d values." % N_VALUES)

            TOTAL_VALUES += (self.NBBBS + self.LIMITR) * 2
            self.log('%d actual values required according to EQDSK G format.' %
                     TOTAL_VALUES)

            self.RBBBS = []
            self.ZBBBS = []
            for i in range(self.NBBBS):
                self.RBBBS.append(P2[_I])
                self.ZBBBS.append(P2[_I + 1])
                _I += 2

            self.RLIM = []
            self.ZLIM = []
            for i in range(self.LIMITR):
                self.RLIM.append(P2[_I])
                self.ZLIM.append(P2[_I + 1])
                _I += 2
        except Exception as e:
            log.error("Failed to read EQDSK")
            self.successfullRead = False
            return False
        self.successfullRead = True
        return True

    def readOneLine(self, line: str) -> list:
        """Reads the float values from one line.

        Returns:
            out (Union[int, List[Float]]): Returns -1 if a line contains an
                alpha at the beginning, else it returns a list of values.
        """
        pattern = r'-?[0-9]{1}\.[0-9]{9}[eE][-+]?[0-9]{2}'

        result = re.findall(pattern, line)

        if line.lstrip()[0].isalpha():
            return -1

        if len(result) == 5:
            # 5 float values on line.
            out = [float(el) for el in result]
            return out
        else:
            # Manually check the values
            out = []
            for el in line.split():
                try:
                    out.append(float(el))
                except Exception as e:
                    pass
            if len(out) == 5:
                # 5 values. All ok
                return out
            else:
                # Must be a comment line
                return -1

    def readValues(self, LINES):
        """Reads values from EQDSK lines.
        """
        # Part one consists of values before the section of limiter and
        # boundary geometry. Basically meaning reading values until finding
        # two integer values.
        part_1 = []
        # Part two are all values regarding plasma and boundary geometry.
        # HOPEFULLY, the lines start with two INTEGERS and not somewhere in
        # between.
        part_2 = []

        # # Pattern for reading integer values separated with spaces
        # integerPattern = r'(?:^|(?<=\s))-?[0-9]+(?:$|(?=\s))'
        # intPatt = re.compile('.*' + integerPattern + '.*')
        # # Pattern for reading scientific numbers
        # scientificPattern = r'[+-]?[0-9]\.[0-9]{6,9}[eE][+-][0-9]{2}'

        # pattern = u'(?:%s|%s)' % (integerPattern, scientificPattern)
        # genPatt = re.compile(pattern)
        SWITCH = 1

        for line in LINES:
            if SWITCH and self.intPatt.match(line):
                SWITCH = 0
            result = self.genPatt.findall(line)
            if SWITCH:
                part_1 += [float(el) for el in result]
            else:
                part_2 += [float(el) for el in result]

        return part_1, part_2

    def readOk(self) -> bool:
        return self.successfullRead

    def generateText(self) -> str:
        """Generates an EQDSK format text.
        """
        # return self.eqdskString

        # The problem with EQDSK G format is that is no longer a consistent
        # standard. Therefore the utility eqdsk.py is mainly used for reading
        # the value, maybe validating the format and using the read value
        # to plot data. Therefore when we want to save the contents of the
        # eqdskString, the eqdsk is no longer generated by the EQDSK G format
        # but the contents of the loaded eqdsk file is just printed.

        #---------------------------------------------------------------------
        if not self.successfullRead:
            return self.eqdskString

        self.generatedText = ''

        self.generatedText += self.HEADER + '\n'

        # First line
        DUMMY = 0.0
        self._PerLineCount = 1


        self.createRealLine(self.RDIM, self.ZDIM, self.RCENTR,
                            self.RLEFT, self.ZMID)

        self.createRealLine(self.RMAXIS,
                            self.ZMAXIS,
                            self.SIMAG, self.SIBRY,
                            self.BCENTR)

        self.createRealLine(self.CURRENT, self.SIMAG, DUMMY,
                            self.RMAXIS, DUMMY)

        self.createRealLine(self.ZMAXIS, DUMMY,
                            self.SIBRY, DUMMY, DUMMY)

        self.createRealLine(*self.FPOL)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.PRES)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.FFPRIM)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.PPRIME)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        for psi_width in self.PSIRZ:
            self.createRealLine(*psi_width)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.QPSI)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'

        self.createIntLine(self.NBBBS, self.LIMITR)
        self.generatedText += '\n'
        self._PerLineCount = 1

        alternating_array = [None] * 2 * self.NBBBS
        alternating_array[::2] = self.RBBBS
        alternating_array[1::2] = self.ZBBBS
        self.createRealLine(*alternating_array)
        self.generatedText += '\n'

        self._PerLineCount = 1
        alternating_array = [None] * 2 * self.LIMITR
        alternating_array[::2] = self.RLIM
        alternating_array[1::2] = self.ZLIM
        self.createRealLine(*alternating_array)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'

        return self.generatedText

    def createRealLine(self, *args) -> None:
        Format = '%.9E'
        for el in args:
            stringReal = Format % el
            self.generatedText += ' ' * (16 - len(stringReal)) + stringReal
            if self._PerLineCount % 5 == 0:
                self.generatedText += '\n'
            self._PerLineCount += 1

    def createIntLine(self, *args):
        for el in args:
            self.generatedText += '   ' + str(el)

    def mmMoveInR(self, mm):
        """ Radially moves the values of eqdsk for *mm* millimeters.

        All R relevant quantities:
            RLIM
            RBBBS
            RLEFT
            RMAXIS

        Note that these quantities are not necessarily in millimeters!
        """

        if not self.successfullRead:
            return

        orderLim = '%.1E' % max(self.RLIM)
        orderLim = int(orderLim[-2:])

        scale = 1.0
        if orderLim == 0:
            # Meters
            scale = 1e-3
        elif orderLim == 2:
            # Centimeters
            scale = 1e-2
        elif orderLim == 3:
            # Millimeters
            scale = 1
        else:
            # print('Strange units as RLIM!')
            return

        for i in range(len(self.RLIM)):
            self.RLIM[i] += scale * mm


        orderBoundary = '%.1E' % max(self.RBBBS)
        orderBoundary = int(orderBoundary[-2:])

        scale = 1.0

        if orderBoundary == 0:
            scale = 1e-3
        elif orderBoundary == 2:
            scale = 1e-2
        elif orderBoundary == 3:
            scale = 1
        else:
            # print('Strange Units at RBBBS!')
            return

        for i in range(len(self.RBBBS)):
            self.RBBBS[i] += scale * mm

        orderLeftR = '%.1E' % self.RLEFT
        orderLeftR = int(orderLeftR[-2:])

        scale = 1.0

        if orderLeftR == 0:
            scale = 1e-3
        elif orderLeftR == 2:
            scale = 1e-2
        elif orderLeftR == 3:
            scale = 1
        else:
            # print('Strange Units at RLEFT!')
            return

        self.RLEFT += scale * mm

        orderMagAxisR = '%.1E' % self.RLEFT
        orderMagAxisR = int(orderMagAxisR[-2:])

        scale = 1.0

        if orderMagAxisR == 0:
            scale = 1e-3
        elif orderMagAxisR == 2:
            scale = 1e-2
        elif orderMagAxisR == 3:
            scale = 1
        else:
            # print('Strange Units at RMAXIS!')
            return

        self.RMAXIS += scale * mm

    def fluxOffset(self, offset):
        """Offsetting the flux quantities.
        SIMAG += offset
        SIBRY = offset
        PSIRZ += offset
        """
        if not self.successfullRead:
            return

        self.SIMAG += offset
        self.SIBRY = offset

        X = len(self.PSIRZ)
        Y = len(self.PSIRZ[0])

        for i in range(X):
            for j in range(Y):
                self.PSIRZ[i][j] += offset

    def convertTo_m(self, a):
        """Determine the size unit. If the units are in mm, divide the
        values by 1000.

        Arguments:
            a (array): Array of floats holding the R or Z values.
        """
        import numpy as np
        order = int(np.log10(max([abs(el) for el in a])))
        if order == 2:
            # Hopefully centimeters
            return [el / 100.0 for el in a]
        elif order == 3:
            # Hopefully mm
            return [el / 1000.0 for el in a]
        else:
            return a

    def getRBBBS(self):
        """Gets the R points for LCFS.

        """
        return self.convertTo_m(self.RBBBS)

    def getZBBBS(self):
        """Gets the Z points for LCFS.

        Determine the size unit. If the units are in meters, multiply the
        values by 1000.
        """
        return self.convertTo_m(self.ZBBBS)

    def getRLIM(self):
        """Gets the R points for wall silhouette.
        """
        return self.convertTo_m(self.RLIM)

    def getZLIM(self):
        """Gets the Z points for wall silhouette.
        """
        return self.convertTo_m(self.ZLIM)

    def getCURRENT(self):
        """Returns plasma current in Ampere
        """
        return self.CURRENT

    def getHEADER(self):
        return self.HEADER

    def getPSIRZ(self):
        return self.PSIRZ

    def getRDIM(self):
        return self.RDIM

    def getZDIM(self):
        return self.ZDIM

    def getNW(self):
        return self.NW

    def getNH(self):
        return self.NH

    def getLIMITR(self):
        return self.LIMITR

    def getNBBBS(self):
        return self.NBBBS

    def getRMAXIS(self):
        return self.RMAXIS

    def getZMAXIS(self):
        return self.ZMAXIS

    def getFPOL(self):
        return self.FPOL

    def getPRES(self):
        return self.PRES

    def getFFPRIM(self):
        return self.FFPRIM

    def getPPRIME(self):
        return self.PPRIME

    def getQPSI(self):
        return self.QPSI

    def getSIMAG(self):
        return self.SIMAG

    def getSIBRY(self):
        return self.SIBRY

    def getRCENTR(self):
        return self.RCENTR

    def getBCENTR(self):
        return self.BCENTR

    def getRLEFT(self):
        return self.RLEFT

    def getZMID(self):
        return self.ZMID

    def setHEADER(self, new):
        self.HEADER = new

    def setPSIRZ(self, new):
        self.PSIRZ = new

    def setCURRENT(self, new):
        self.CURRENT = new

    def setRDIM(self, new):
        self.RDIM = new

    def setZDIM(self, new):
        self.ZDIM = new

    def setNW(self, new):
        self.NW = new

    def setNH(self, new):
        self.NH = new

    def setLIMITR(self, new):
        self.LIMITR = new

    def setRLIM(self, new):
        self.RLIM = new

    def setZLIM(self, new):
        self.ZLIM = new

    def setNBBBS(self, new):
        self.NBBBS = new

    def setRBBBS(self, new):
        self.RBBBS = new

    def setZBBBS(self, new):
        self.ZBBBS = new

    def setRMAXIS(self, new):
        self.RMAXIS = new

    def setZMAXIS(self, new):
        self.ZMAXIS = new

    def setFPOL(self, new):
        self.FPOL = new

    def setPRES(self, new):
        self.PRES = new

    def setFFPRIM(self, new):
        self.FFPRIM = new

    def setPPRIME(self, new):
        self.PPRIME = new

    def setQPSI(self, new):
        self.QPSI = new

    def setSIMAG(self, new):
        self.SIMAG = new

    def setSIBRY(self, new):
        self.SIBRY = new

    def setRCENTR(self, new):
        self.RCENTR = new

    def setBCENTR(self, new):
        self.BCENTR = new

    def setRLEFT(self, new):
        self.RLEFT = new

    def setZMID(self, new):
        self.ZMID = new


def createEqdskFromSlice(slice, vacuum_toroidal_field, HEADER="") -> EQDSKIO:
    """From a slice of the equilibrium IDS, create an EQDSKIO object.
    """
    import numpy as np
    eqObj = EQDSKIO()

    profile2d = slice.profiles_2d[0]

    # Rectangular grid
    R = profile2d.grid.dim1
    Z = profile2d.grid.dim2

    Nw = R.shape[0]
    Nh = Z.shape[0]

    eqObj.setNW(Nw)
    eqObj.setNH(Nh)

    # Setting R,Z grid
    eqObj.setRLEFT(R[0])
    eqObj.setRDIM(R[-1] - R[0])
    eqObj.setZMID((Z[0] + Z[-1]) / 2)
    eqObj.setZDIM(Z[-1] - Z[0])
    psi = profile2d.psi / (2*np.pi)
    eqObj.setPSIRZ(psi.T)

    profile1d = slice.profiles_1d

    def densify(f, N):
        """Extends the array to N if shape is different.
        """
        from scipy.interpolate import interp1d
        if f.shape[0] == 0:
            return np.array([0 for i in range(N)])

        if f.shape[0] == N:
            return f

        x = np.linspace(0, 1, f.shape[0])
        i = interp1d(x, f)
        newX = np.linspace(0, 1, N)
        newF = i(newX)
        return newF

    # print(profile1d.f.shape, profile1d.pressure.shape,
    #       profile1d.dpressure_dpsi.shape)
    FPOL = densify(profile1d.f,Nw)

    eqObj.setFPOL(FPOL)
    eqObj.setPRES(densify(profile1d.pressure, Nw))
    eqObj.setPPRIME(densify(profile1d.dpressure_dpsi, Nw) / (2 * np.pi))
    eqObj.setFFPRIM(densify(profile1d.f_df_dpsi, Nw) / (2 * np.pi))
    eqObj.setQPSI(densify(profile1d.q, Nw)) # empty

    globQuant = slice.global_quantities


    eqObj.setRMAXIS(globQuant.magnetic_axis.r)
    eqObj.setZMAXIS(globQuant.magnetic_axis.z)
    eqObj.setSIBRY(globQuant.psi_boundary / (2 * np.pi))
    eqObj.setSIMAG(globQuant.psi_axis / (2 * np.pi))
    eqObj.setCURRENT(globQuant.ip)
    eqObj.setRCENTR(globQuant.magnetic_axis.r)
    BCENTR = vacuum_toroidal_field.b0[0]
    # print(f"B_t at centre: {BCENTR} T")
    current = np.abs(globQuant.ip)
    current = current / 1e6
    current = round(current)
    # print(f"Drsep: {slice.boundary_separatrix.gap[-1].value}")
    if len(slice.boundary_separatrix.gap):
        drsep = slice.boundary_separatrix.gap[-1].value
        HEADER = f"{current}MA {HEADER} dsep={100 * drsep:.02f}cm - 3 {Nw} {Nh}"
    else:
        HEADER = f"{current}MA {HEADER} - 3 {Nw} {Nh}"
    eqObj.setHEADER(HEADER)
    eqObj.setBCENTR(BCENTR)

    boundary = slice.boundary
    eqObj.setRBBBS(boundary.outline.r)
    eqObj.setZBBBS(boundary.outline.z)
    eqObj.setNBBBS(len(boundary.outline.r))


    return eqObj

def getBackUpIMASWallIds(shot=116000, run=4, db_name="ITER_MD",
        user_name="public"):
    """Gets the default Wall IDS machine description from the ITER ITER_MD
    machine description database. 116000/4
    """

    import imas
    import imas.imasdef

    db_entry = imas.DBEntry(shot=shot, run=run, db_name=db_name,
            user_name=user_name, backend_id=imas.imasdef.MDSPLUS_BACKEND)
    db_entry.open()

    return db_entry.get("wall")

def addWallDescriptionToEqdsk(eqObj, idsWall):
    """From the wall IDS take the wall silhouette points and save it to the
    EQDSKIO object.
    """
    if len(idsWall.description_2d) == 0:
        idsWall = getBackUpIMASWallIds()

    description2d = idsWall.description_2d[0]
    eqObj.setRLIM(description2d.limiter.unit[0].outline.r)
    eqObj.setZLIM(description2d.limiter.unit[0].outline.z)
    eqObj.setLIMITR(len(description2d.limiter.unit[0].outline.r))

description="""Create an EQDSK G file from an ids.equilibrium.time_slice.

The name of the generated file:

shot_$_run_$_t_$s.eqdsk

The file is outputted wherever the command is run.

"""

parser = argparse.ArgumentParser(description='Create EQDSK from IDS time slice')

parser.add_argument('-u', '--user',
                    metavar='user',
                    type=str,
                    help='argument for user. Default public',
                    default='public')
parser.add_argument('-d', '--device',
                    metavar='device',
                    type=str,
                    help='argument for device. Default ITER',
                    default="ITER")
parser.add_argument('-v', '--version',
                    metavar='version',
                    type=str,
                    default='3',
                    help='argument for version. Default 3')
parser.add_argument('-s', '--shot',
                    type=int,
                    metavar='shot',
                    help='argument for shot. Required',
                    required=True)
parser.add_argument('-r', '--run',
                    type=int,
                    metavar='run',
                    help='argument for run. Required',
                        required=True)
parser.add_argument('-ts', '--time_start', default=0.0,
                    type=float, metavar='time start',
                    help='Start of time interval.')
parser.add_argument('-te', '--time_end', default=float("inf"),
                    type=float, metavar='time end',
                    help='End of time interval.')
parser.add_argument('-b', '--backend', type=str, metavar="backend type",
                default="MDSPLUS")
args = parser.parse_args()

user = args.user
device = args.device
shot = args.shot
run = args.run
version = args.version
# t_s = 399.927598 # time slice
DINA_SCENARIO_NAME="DT-DINA2020-04"

if args.backend == "MDSPLUS":
    backend_id = imas.imasdef.MDSPLUS_BACKEND
elif args.backend == "HDF5":
    backend_id = imas.imasdef.HDF5_BACKEND
else:
    raise Exception(f"Wrong -b or --backend option: '{args.backend}'")

entry = imas.DBEntry(shot=shot, run=run, user_name=user,
                     backend_id=backend_id,
                     db_name=device, data_version=args.version)
entry.open()
# imas_obj = imas.ids(shot, run)
# imas_obj.open_env(user, device, version)

# ids_equilibrium = imas_obj.equilibrium
# ids_equilibrium = entry.get("equilibrium")

# time_slices = args.time_slices[:]

# ids_wall = imas_obj.wall
# ids_wall.get()
# ids_wall = entry.get("wall")
ids_wall = getBackUpIMASWallIds()


# Check if time is homogeneous
equil_properties = entry.partial_get("equilibrium", "ids_properties")
if equil_properties.homogeneous_time:
    # Get the times from here
    ids_times = entry.partial_get("equilibrium", "time")
else:
    # Get it from summary
    summary = entry.get("summary")
    ids_times = summary.time


for ts in ids_times:
    if ts < args.time_start or ts > args.time_end:
        continue
    print(f"Processing {ts}")
    # ids_equilibrium.getSlice(ts, imas.imasdef.CLOSEST_SAMPLE)
    ids_equilibrium = entry.get_slice("equilibrium", ts, imas.imasdef.CLOSEST_SAMPLE)

    # Actual time depends on what is the closest sample, do NOT trust linear
    # interpolation!!


    slice = ids_equilibrium.time_slice[0]

    HEADER_COMMENT=f"IDS shot={shot} run={run} t={ts:.3f}s"
    eqObj = createEqdskFromSlice(slice, ids_equilibrium.vacuum_toroidal_field,
                                           HEADER_COMMENT)
    # ids_wall = imas_obj.wall
    # ids_wall.get()
    addWallDescriptionToEqdsk(eqObj, ids_wall)

    # SuccessfullRead is needed to set to 1 as we manually populate the data
    # and not read from an existing EQDSK G file. Without this, the
    # generateText function would not return anything.
    eqObj.successfullRead = True
    text = eqObj.generateText()
    f = open(f'shot_{shot}_run_{run}_t_{ts:.3f}s.eqdsk', 'w')
    f.write(text)
    f.close()
