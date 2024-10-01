import datetime
from pathlib import Path
import FreeCAD as App # type: ignore


def readInput(self):
    
    paraFile = self.form[0].setupFile.text()

    # Read parameters from file
    with open(Path(paraFile), "r") as f:
        for line in f:
            if "radial_growth_rule" in line:
                radial_growth_rule = line.split("=")[1].strip()
            elif "r_min" in line:
                r_min = line.split("=")[1].strip()
            elif "r_max" in line:
                r_max = line.split("=")[1].strip()
            elif "nrings" in line:
                nrings = line.split("=")[1].strip()

            elif "cellsize_early" in line:
                cellsize_early = line.split("=")[1].strip()
            elif "cellsize_late" in line:
                cellsize_late = line.split("=")[1].strip()
            elif "cellwallthickness_early" in line:
                cellwallthickness_early = line.split("=")[1].strip()
            elif "cellwallthickness_late" in line:
                cellwallthickness_late = line.split("=")[1].strip()

            elif "skeleton_density" in line:
                skeleton_density = line.split("=")[1].strip()
            elif "random_noise" in line:
                random_noise = line.split("=")[1].strip()

            # elif "iter_max" in line:
            #     iter_max = line.split("=")[1].strip()

            elif "box_center" in line:
                box_center = line.split("=")[1].strip()
            elif "box_height" in line:
                box_height = line.split("=")[1].strip()
            elif "box_width" in line:
                box_width = line.split("=")[1].strip()
            elif "box_depth" in line:
                box_depth = line.split("=")[1].strip()

            elif "nsegments" in line:
                nsegments = line.split("=")[1].strip()
            elif "theta_min" in line:
                theta_min = line.split("=")[1].strip()
            elif "long_connector_ratio" in line:
                long_connector_ratio = line.split("=")[1].strip()

            elif "Uinf" in line:
                Uinf = line.split("=")[1].strip()
            elif "a1" in line:
                a1 = line.split("=")[1].strip()
            elif "a2" in line:
                a2 = line.split("=")[1].strip()
            elif "m1" in line:
                m1 = line.split("=")[1].strip()
            elif "m2" in line:
                m2 = line.split("=")[1].strip()

            elif "x_notch_size" in line:
                x_notch_size = line.split("=")[1].strip()
            elif "y_notch_size" in line:
                y_notch_size = line.split("=")[1].strip()

            elif "precrackFlag" in line:
                precrackFlag = line.split("=")[1].strip()
                if precrackFlag in ['on','On','Y','y','Yes','yes'] and "precrack_size" in line:
                    precrack_size = line.split("=")[1].strip()

            elif "boundaryFlag" in line:
                boundaryFlag = line.split("=")[1].strip()
            elif "box_shape" in line:
                box_shape = line.split("=")[1].strip()
                
            elif "merge_operation" in line:
                merge_operation = line.split("=")[1].strip()
                if merge_operation in ['on','On','Y','y','Yes','yes'] and "merge_tol" in line:
                    merge_tol = line.split("=")[1].strip()
                
            elif "stlFlag" in line:
                stlFlag = line.split("=")[1].strip()
            elif "inpFlag" in line:
                inpFlag = line.split("=")[1].strip()
            elif "inpType" in line:
                inpType = line.split("=")[1].strip()

    # Write parameters to input panel
        self.form[0].radial_growth_rule.setCurrentText(radial_growth_rule)
        self.form[0].r_min.setText(r_min)
        self.form[0].r_max.setText(r_max)
        self.form[0].nrings.setText(nrings)

        self.form[0].cellsize_early.setText(cellsize_early)
        self.form[0].cellsize_late.setText(cellsize_late)
        self.form[0].cellwallthickness_early.setText(cellwallthickness_early)
        self.form[0].cellwallthickness_late.setText(cellwallthickness_late)

        self.form[0].skeleton_density.setText(skeleton_density)
        self.form[0].random_noise.setText(random_noise)

        # self.form[0].iter_max.setText(iter_max)

        self.form[0].box_center.setText(box_center)
        self.form[0].box_height.setText(box_height)
        self.form[0].box_width.setText(box_width)
        self.form[0].box_depth.setText(box_depth)
        
        self.form[0].nsegments.setText(nsegments)
        self.form[0].theta_min.setText(theta_min)
        self.form[0].long_connector_ratio.setText(long_connector_ratio)

        self.form[0].x_indent_size.setText(x_notch_size)   
        self.form[0].y_indent_size.setText(y_notch_size)

        self.form[0].precrackFlag.setCurrentText(precrackFlag)
        if precrackFlag in ['on','On','Y','y','Yes','yes']:
            self.form[0].precrack_size.setText(precrack_size)

        self.form[0].boundaryFlag.setCurrentText(boundaryFlag)
        self.form[0].box_shape.setCurrentText(box_shape)
        
        self.form[0].merge_operation.setCurrentText(merge_operation)
        # if merge_operation in ['on','On','Y','y','Yes','yes']:
        #     self.form[0].merge_tol.setText(merge_tol)

        self.form[0].stlFlag.setCurrentText(stlFlag)
        self.form[0].inpFlag.setCurrentText(inpFlag)
        self.form[0].inpType.setCurrentText(inpType)