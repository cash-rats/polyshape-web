import { DownloadIcon } from '@radix-ui/react-icons';
import * as VisuallyHidden from "@radix-ui/react-visually-hidden";
import {
  Dialog,
  DialogContent,
  DialogTrigger,
  DialogTitle,
} from "~/components/ui/dialog";
import { useTrianglifyPreview } from '~/routes/hooks/use-trianglify-preview';
import { Dimensions } from "~/types";

interface ExportDialogProps {
  triangifyPattern: any;
  dimensions: Dimensions;
}

export default function ExportDialog({ triangifyPattern, dimensions }: ExportDialogProps) {
  console.log('ExportDialog', triangifyPattern);
  const { previewRef } = useTrianglifyPreview(triangifyPattern, dimensions);

  return (
    <Dialog>
      <DialogTrigger asChild>
        <button
          type="button"
          className="flex items-center px-4 py-2 bg-white text-[#578E7E] rounded-md hover:bg-gray-100"
        >
          <DownloadIcon className="w-5 h-5 stroke-current text-[#AEEA94]" />
          <span className="ml-2">Export</span>
        </button>
      </DialogTrigger>
      <DialogContent aria-describedby={undefined} className="max-w-[600px] w-full p-0">
        <VisuallyHidden.Root>
          <DialogTitle>Export Pattern</DialogTitle>
        </VisuallyHidden.Root>

        {/* Pattern preview*/}
        <div
          ref={previewRef}
          className="w-full h-full flex items-center justify-center bg-[#e8e8e8] p-7"
        />

        <div className="grid gap-4 py-4">
          hello
        </div>
      </DialogContent>
    </Dialog>
  );
}
