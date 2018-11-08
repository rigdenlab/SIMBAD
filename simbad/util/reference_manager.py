"""Code to manage references"""

__author__ = "Adam Simpkin, Jens Thomas & Felix Simkovic"
__date__ = "7 November 2018"

from collections import OrderedDict
import copy
from enum import Enum
import os


class ReferenceManager():

    #Section Names
    class SECTIONS(Enum):
        __order__ = 'GENERAL ROTATION MR REFINEMENT'
        GENERAL = 'General'
        ROTATION = 'Rotation Search'
        MR = 'Molecular Replacement'
        REFINEMENT = 'Refinement'

    SEC_TAG = 'h3'

    def __init__(self, d):
        self.references = {}
        self.ordered_labels = []
        self.citation_file_path = None
        self.section_labels = OrderedDict()
        self.setup_references()
        self.setup_sections(d)

    def setup_references(self):
        ref_fname = os.path.join(os.environ['CCP4'], "share", "simbad", "static", "simbad.bib")
        if not os.path.isfile(ref_fname):
            msg = "Cannot find BibTex file containing references. " \
                  "Please determine them yourself and cite AMPLE."
            return msg
        article = {}
        entry = False
        with open(ref_fname, "r") as fhin:
            for line in fhin.readlines():
                line = line.strip()
                if not line:
                    continue
                elif line.startswith("@"):
                    # Beginning of all BibTex entry blocks
                    entry = True
                    unique_id = line.replace("@article{", "").replace(",", "")
                    article = {'unique_id': unique_id}  # Reset the article dictionary
                elif line == "}":
                    # End of all BibTex entry blocks
                    entry = False
                    self.references[article['label']] = article
                elif entry:
                    # BibTex entry block
                    # Some dirty line handling.
                    # Not very bulletproof but should do for now
                    line = line.replace("{", "").replace("}", "")
                    key, value = [l.strip() for l in line.split("=")]
                    value = value.rstrip(",").replace("\"", "")
                    # Save the data to the article entry
                    article[key] = value
        return

    def setup_sections(self, d):
        # Create default lists
        for section in self.SECTIONS:
            self.section_labels[section] = []
            # Build up list of program reference labels, ordered by sections
            for section in self.SECTIONS:
                if section == self.SECTIONS.GENERAL:
                    self.section_labels[section] = ['SIMBAD', 'CCP4', 'CCTBX']
                elif section == self.SECTIONS.ROTATION:
                    labels = []
                    if d.get('-rot_program') == 'phaser':
                        labels.append('PHASER_ROT')
                    else:
                        labels.append('AMORE')
                    self.section_labels[section] = labels
                elif section == self.SECTIONS.MR:
                    labels = []
                    if d.get('-mr_program') == 'phaser':
                        labels.append('PHASER')
                    else:
                        labels.append('MOLREP')
                    self.section_labels[section] = labels
                elif section == self.SECTIONS.REFINEMENT:
                    self.section_labels[self.SECTIONS.REFINEMENT] = ['REFMAC']

        # Generate ordered list of all relevant reference labels
        for section in self.SECTIONS:
            if section in self.section_labels:
                self.ordered_labels += self.section_labels[section]
        return

    @property
    def methods_as_html(self):
        html = "<p>This section lists the programs and algorithms that were used in this job and the references that should be cited. " + \
               "Numbers in superscript next to program/reference names refer to the number of the program reference in the overall list of references.</p>"
        for section in self.SECTIONS:
            if section == self.SECTIONS.GENERAL:
                html += '<p>The first 2 references should be cited in all cases.</p>'
            elif section == self.SECTIONS.ROTATION and len(self.section_labels[self.SECTIONS.ROTATION]):
                standfirst = '<p>The rotation searches were carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.ROTATION, standfirst)
            elif section == self.SECTIONS.MR and len(self.section_labels[self.SECTIONS.MR]):
                standfirst = '<p>Molecular Replacement was carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.MR, standfirst)
            elif section == self.SECTIONS.REFINEMENT and len(self.section_labels[self.SECTIONS.REFINEMENT]):
                standfirst = '<p>Refinement of the MR solutions carried out with the following programs:</p>'
                html += self._methods_section_html(self.SECTIONS.REFINEMENT, standfirst)
        return html

    def _methods_section_html(self, section, standfirst):
        mysec = self.SECTIONS(section)
        html = '<{}>{}</{}>'.format(self.SEC_TAG, mysec.value, self.SEC_TAG)
        html += standfirst
        html += '<ul>'
        for label in self.section_labels[mysec]:
            html += "<li>{}<sup>{}</sup></li>".format(label, self.ordered_labels.index(label) + 1)
        html += "</ul>"
        return html

    @property
    def citations_as_html(self):
        html = '<{}>References</{}>'.format(self.SEC_TAG, self.SEC_TAG)
        html += '<ol>'
        template_txt = "<li> {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]</li>"
        for label in self.ordered_labels:
            ref = copy.copy(self.references[label])
            ref['author'] = ref['author'].split(" and ")[0].split(",")[0] + " et al."
            ref['pages'] = ref['pages'].replace("--", "-")
            html += template_txt.format(**ref)
        html += '</ol>'
        return html

    @property
    def citations_as_text(self):
        txt = """A number of programs and algorithms were used within the this run of AMPLE.
    The following is a list of citations for this run:
    {0}
    """.format(self.citation_list_as_text)
        if self.citation_file_path:
            txt += """
    A bibtex file with these references has been saved to the following file:
    {0}
    """.format(self.citation_file_path)
            return txt

    @property
    def citation_list_as_text(self):
        template_txt = "* {author} ({year}). {title}. {journal} {volume}({number}), {pages}. [doi:{doi}]"
        text = ""
        for label in self.ordered_labels:
            ref = copy.copy(self.references[label])
            ref['author'] = ref['author'].split(" and ")[0].split(",")[0] + " et al."
            ref['pages'] = ref['pages'].replace("--", "-")
            text += template_txt.format(**ref) + os.linesep * 2
        return text

    def save_citations_to_file(self, work_dir):
        # =========================================================================
        # Somewhat a template of how we want to write each article in BibTex format
        # =========================================================================
        template_bib = "@article{{{unique_id},{sep}author = {{{author}}},{sep}doi = {{{doi}}},{sep}" \
                       "journal = {{{journal}}},{sep}number = {{{number}}},{sep}pages = {{{pages}}},{sep}" \
                       "title = {{{{{title}}}}},{sep}volume = {{{volume}}},{sep}year = {{{year}}},{sep}}}{sep}"
        references_bib = [template_bib.format(sep=os.linesep, **self.references[l]) for l in self.ordered_labels]
        ref_fname = os.path.join(work_dir, "references.bib")
        with open(ref_fname, "w") as fhout:
            fhout.write(os.linesep.join(references_bib))
        self.citation_file_path = ref_fname
        return ref_fname