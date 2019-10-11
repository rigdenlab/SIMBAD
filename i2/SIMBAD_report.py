"""
        SIMBAD_report.py: CCP4 GUI Project

        This library is free software: you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public License
        version 3, modified in accordance with the provisions of the
        license to address the requirements of UK law.

        You should have received a copy of the modified GNU Lesser General
        Public License along with this library.  If not, copies may be
        downloaded from http://www.ccp4.ac.uk/ccp4license.php

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU Lesser General Public License for more details.
        """


import os
import re

if __name__ == "__main__":
    import sys

    ccp4 = os.environ["CCP4"]
    sys.path.append(os.path.join(ccp4, "share", "ccp4i2", "report"))
    sys.path.append(os.path.join(ccp4, "share", "ccp4i2", "core"))
    sys.path.append(os.path.join(ccp4, "lib", "python2.7", "site-packages"))

from lxml import etree as ET
from report.CCP4ReportParser import Report

from simbad.util import SIMBAD_DIRNAME, SIMBAD_PYRVAPI_SHAREDIR


class SIMBAD_report(Report):
    TASKNAME = "SIMBAD"
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        repdir = os.path.join(jobInfo.get("fileroot", None), SIMBAD_DIRNAME, SIMBAD_PYRVAPI_SHAREDIR)
        self.get_tables_as_elements(repdir)
        # print("JMHT WRITING REPORT %s" % self.e1_dict)
        self.addDiv(style="clear:both;")
        for e1 in xmlnode:
            # Process each tab separately
            if e1.tag == "tab":
                self.report_section(e1, self)
        return

    def get_tables_as_elements(self, repdir):
        """Get tables as xmltree elements by parsing task.tsk file and .table files"""
        try:
            t1_list = list()
            with open(os.path.join(repdir, "task.tsk")) as istream:
                # print("JMHT CHECKING task.tsk %s\n" % os.path.join(repdir, 'task.tsk'))
                for s1 in re.findall("<table .+?</table>", istream.read(), re.S):
                    t1 = ET.fromstring(s1)
                    if len(t1):
                        t1_list.append(t1)
            for f1 in os.listdir(repdir):
                if f1.endswith(".table"):
                    t1 = ET.parse(os.path.join(repdir, f1)).getroot()
                    if len(t1):
                        t1_list.append(t1)
            self.e1_dict = dict()
            for t1 in t1_list:
                tid = t1.get("id", None)
                if tid and tid.endswith("-grid"):
                    tags = [t2.tag for t2 in t1]
                    if tags == ["thead", "tbody"]:
                        assert len(t1) == 2
                        e1 = t1
                    else:
                        tset = set(tags)
                        tag = tset.pop()
                        assert not tset and tag == "tr"
                        e1 = ET.Element("table")
                        e1.append(t1)
                        e1.attrib.update(t1.attrib)
                        t1.attrib.clear()
                        t1.tag = "tbody"
                    for e2 in e1.iter():
                        e2.attrib.pop("class", None)
                    e1.find("tbody").set("class", "fancy")
                    self.e1_dict[tid[:-5]] = e1
            if len(self.e1_dict.keys()):
                return True
            return False
        except Exception as e:
            print "EXCEPTION: {0}".format(e)
            return

    def report_section(self, e1, r0, sort=False):
        """

        """
        elems = list()
        title = "Untitled"
        state = False
        cou = 0
        # print("Processing tag %s id %s\n%s" % (e1.tag, e1.get('id'),ET.tostring(e1)))
        for e2 in e1:
            row = e2.get("row", "_")
            col = e2.get("col", "_")
            if row.isdigit():
                row = int(row)
            if col.isdigit():
                col = int(col)
            if e2.get("id") or e2.tag == "text":
                elems.append([row, col, e2])
                if e2.tag == "table":
                    cou += 1
            elif e2.tag == "name":
                title = e2.text.strip()
            elif e2.tag == "open":
                state = e2.text.strip() == "true"
        if elems:
            # strip out anything we can't deal with here
            if any([x in title.lower() for x in ["downloads", "log files", "graph"]]):
                return
            # print "GOT ELEMS ",[g[2].get('id') for g in elems],title
            r1 = r0.addFold(label=title, initiallyOpen=state)
            # for row, col, e2 in sorted(grid):
            if sorted:
                elems = sorted(elems)
            for _, _, e2 in elems:
                id2 = e2.get("id")
                # print "PROCESSING ",id2, e2.tag
                if e2.tag == "section":
                    self.report_section(e2, r1)
                elif e2.tag == "table":
                    if id2 and id2 in self.e1_dict:
                        if id2 == "mrbump_table":
                            r1.append("The table below details the Molecular Replacement results from MrBUMP")
                        if cou > 1:
                            r1.append(e2.findtext("legend").strip())
                        r1.append(ET.tostring(self.e1_dict[id2]))
                elif e2.tag == "text":
                    for t in e2.itertext():
                        r1.append(t)
                else:
                    pass


if __name__ == "__main__":
    # Run with no arguments in the CCP4 job directory (the one that holds the SIMBAD directory)
    def test2():
        import argparse

        parser = argparse.ArgumentParser(
            description="test of morda report generator", formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            "-w",
            "--wrkdir",
            help="""a directory, containing the subdirectory
                report/ generated by rvapi""",
            default=".",
            metavar="<dir>",
        )
        parser.add_argument(
            "-i", "--xml", help="input xml-file generated previously by rvapi", default="program.xml", metavar="<file>"
        )
        parser.add_argument(
            "-o", "--html", help="output html-file, a report file for i2", default="areport.html", metavar="<file>"
        )
        opt = parser.parse_args()
        xmlnode = ET.parse(opt.xml).getroot()
        jobInfo = dict(fileroot=os.path.abspath(opt.wrkdir))
        report = SIMBAD_report(xmlnode, jobInfo)
        if len(report.errReport):
            print "ERROR REPORT"
            print report.errReport.report()

        htmlbase = "file://" + os.environ["CCP4"] + "/share/ccp4i2/docs/report_files"
        htmlstr = ET.tostring(report.as_etree(htmlBase=htmlbase))
        with open(opt.html, "w") as ostream:
            print >> ostream, htmlstr.replace("><", ">\n<")

    test2()


# #from CCP4ReportParser import Report
# # class AMPLE_report(Report):
# #      # Specify which gui task and/or pluginscript this applies to
# #      TASKNAME = 'AMPLE'
# #      RUNNING = False
# #      def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
# #              Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
# #              clearingDiv = self.addDiv(style="clear:both;")
# #              self.addDefaultReport(self)
# #
# #      def addDefaultReport(self, parent=None):
# #              if parent is None: parent=self
# #              if len(self.xmlnode.xpath("LogText")) > 0:
# #                      newFold = parent.addFold(label="Log text", initiallyOpen=True)
# #                      newFold.addPre(text = self.xmlnode.xpath("LogText")[0].text)
