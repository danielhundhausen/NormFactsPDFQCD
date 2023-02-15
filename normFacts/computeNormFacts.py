#!/nfs/dust/cms/user/hundhad/anaconda3/envs/py39/bin/python
from collections import defaultdict
import glob
import json
import os

import uproot


SAMPLES = [""]


class NormComputer:

    QCD_SCALES = [
        "h_murmuf_upup",
        "h_murmuf_upnone",
        "h_murmuf_noneup",
        "h_murmuf_downdown",
        "h_murmuf_downnone",
        "h_murmuf_nonedown"
    ]

    QCD_NAMES = [
        "murmuf_upup",
        "mur_upnone",
        "muf_noneup",
        "murmuf_downdown",
        "mur_downnone",
        "muf_nonedown"
    ]

    def __init__(self):
        self.preselection_output_path = "."

    def read_histos(self, files_list):
        _nominal = 0.0
        _qcd_variations = defaultdict(list)
        _pdf_variations = defaultdict(list)

        for _file in files_list:
            # Skip HDAMP
            if (('hdamp' in _file)):
                continue
            # Keep TTZ and TTW
            if ((('TuneCP' in _file) and ('TTZ' not in _file))) and ((('TuneCP' in _file) and ('TTW' not in _file))):
                continue
            # Skip TuneCP variations in TT
            if (('TuneCP' in _file) and ('TTTo' in _file)):
                continue
            try:
                with uproot.open(_file) as f:
                    _nom = f['UncNorms/sum_event_weights'].to_numpy()
                    _nominal += _nom[0][0]

                    for _scale in self.QCD_SCALES:
                        _var = f[f'UncNorms/{_scale}'].to_numpy()[0][0]
                        _qcd_variations[_scale].append(_var)

                    for i in range(100):
                        _var = f[f'UncNorms/h_pdf_{i}'].to_numpy()[0][0]
                        _pdf_variations[i].append(_var)
            except Exception as e:
                print(e)
                print(f'Skip {_file}, problematic. Check!!!')

        return _nominal, _qcd_variations, _pdf_variations

    def variations_from_samples(self, sample_specifier):
        nominal = 0.0
        qcd_variations = defaultdict(list)
        pdf_variations = defaultdict(list)

        slist = glob.glob(self.preselection_output_path + "/uhh2.AnalysisModuleRunner.MC." + sample_specifier)
        print(f'------ {sample_specifier} ({len(slist)} files) -----')

        _nominal, _qcd_variations, _pdf_variations = self.read_histos(slist)

        nominal += _nominal

        for _scale in self.QCD_SCALES:
            qcd_variations[_scale].append(sum(_qcd_variations[_scale]))

        for i in range(100):
            pdf_variations[i].append(sum(_pdf_variations[i]))

        return nominal, qcd_variations, pdf_variations

    def _compute_pdf_norm_factors(self, nominal, pdf_variations):
        for key in pdf_variations:
            pdf_variations[key] = sum(pdf_variations[key])

        norm_fact_pdf = []
        for i in range(100):
            norm_fact_pdf.append(nominal / pdf_variations[i])
        return norm_fact_pdf

    def _compute_qcd_norm_factors(self, nominal, qcd_variations):
        for _scale in self.QCD_SCALES:
            qcd_variations[_scale] = sum(qcd_variations[_scale])

        norm_fact_qcd = {}
        for _scale, _w_name in zip(self.QCD_SCALES, self.QCD_NAMES):
            norm_fact_qcd[_w_name] = nominal / qcd_variations[_scale]
        return norm_fact_qcd

    def compute_norms(self):
        """
        Compute the qcd and pdf norms for all samples
        as given by the config module and save to file.
        """
        for sample_specifier in ["DY*", "TT*"]:

            nominal, qcd_variations, pdf_variations = self.variations_from_samples(sample_specifier)

            norm_facts = {
                'PDF': self._compute_pdf_norm_factors(nominal, pdf_variations),
                'QCD': self._compute_qcd_norm_factors(nominal, qcd_variations)
            }

            os.makedirs("norm_facts", exist_ok=True)
            with open(f"norm_facts/{sample_specifier}_norm.json", "w") as f:
                f.write(json.dumps(norm_facts))


if __name__ == "__main__":
    norm_computer = NormComputer()
    norm_computer.compute_norms()

