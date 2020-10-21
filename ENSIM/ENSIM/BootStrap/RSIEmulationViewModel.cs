// ---------------------------------------------------------------------------------------------
// RSIEmulationViewModel.cs
// ---------------------------------------------------------------------------------------------
// Summary :
// ---------------------------------------------------------------------------------------------
// Frédéric POINDRON - 29/06/2017, 10:04
// ---------------------------------------------------------------------------------------------
namespace RSIStandardEmulationApp1.BootStrap
{
    using RSI.Common.Core.ViewModels;
    public class RSIEmulationViewModel : EmulationViewModel
    {
        public override string EmulationName => "RSI";

        protected override bool CanTranslate => false;

        protected override string TranslatorProgramFile => null;
    }
}