using RSI.IndissLike.Controls.Helpers;
using System;
using System.Windows;

namespace RSIStandardEmulationApp1.BootStrap
{
    /// <summary>
    /// Interaction logic for App.xaml
    /// </summary>
    public partial class App
    {
        /// <summary>
        /// Gets the name of the application.
        /// </summary>
        /// <returns></returns>
        protected override string GetApplicationName()
        {
            return "RSIStandardEmulationApp1";
        }

        /// <summary>
        /// Gets the user setting filename.
        /// </summary>
        /// <returns></returns>
        protected override string GetUserSettingFilename()
        {
            return "RSIStandardEmulationApp1.xml";
        }
    }
}
