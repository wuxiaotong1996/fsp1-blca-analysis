@echo off
rem qc_rerun_all.bat — final QC: rerun all supplementary analysis scripts end-to-end
cd /d "D:\HuaweiMoveData\Users\86159\Desktop\a...修改\fsp1_blca_paper"
set RSCRIPT=C:\Users\86159\AppData\Local\Programs\R\R-4.6.1\bin\Rscript.exe
set LOG=results\QC_rerun.log
echo ==== QC rerun START %date% %time% ==== > %LOG%
%RSCRIPT% scripts\analysis1b_downsampling.R >> %LOG% 2>>&1 || echo [FAILED] analysis1b >> %LOG%
%RSCRIPT% scripts\analysis2_estimate_ppcor.R >> %LOG% 2>>&1 || echo [FAILED] analysis2 >> %LOG%
%RSCRIPT% scripts\analysis2b_purity_24class.R >> %LOG% 2>>&1 || echo [FAILED] analysis2b >> %LOG%
%RSCRIPT% scripts\analysis4_ssgsea_24class.R >> %LOG% 2>>&1 || echo [FAILED] analysis4 >> %LOG%
%RSCRIPT% scripts\analysis35_scRNA_lineages.R >> %LOG% 2>>&1 || echo [FAILED] analysis35 >> %LOG%
%RSCRIPT% scripts\analysis6_gse149652_tcell.R >> %LOG% 2>>&1 || echo [FAILED] analysis6 >> %LOG%
%RSCRIPT% scripts\analysis6b_gse135337_aifm2.R >> %LOG% 2>>&1 || echo [FAILED] analysis6b >> %LOG%
echo ==== QC rerun END %date% %time% ==== >> %LOG%
